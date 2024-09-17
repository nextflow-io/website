import sanityClient from '@sanity/client';
import fs from 'fs';
import path from 'path';
import axios from 'axios';
import * as cheerio from 'cheerio';
import { customAlphabet } from 'nanoid'

const nanoid = customAlphabet('0123456789abcdef', 12)

const client = sanityClient({
  projectId: 'o2y1bt2g',
  dataset: 'seqera',
  token: process.env.SANITY_TOKEN,
  useCdn: false,
});

const postsFile = path.join(process.cwd(), 'export.json');
const posts = JSON.parse(fs.readFileSync(postsFile, 'utf8'));

async function downloadImage(url) {
  const response = await axios.get(url, { responseType: 'arraybuffer' });
  return Buffer.from(response.data, 'binary');
}

async function uploadImageToSanity(imageBuffer, filename) {
  return client.assets.upload('image', imageBuffer, { filename });
}

async function replaceImageUrls(content, imageMap) {
  const $ = cheerio.load(content);
  $('img').each((i, elem) => {
    const src = $(elem).attr('src');
    if (src && imageMap[src]) {
      $(elem).attr('src', imageMap[src]);
    }
  });
  return $.html();
}

async function migratePosts() {
  const p = [posts[0], posts[1]]
  for (const post of p) {
    const imageMap = {};
    for (const imageUrl of post.images) {
      try {
        const imageBuffer = await downloadImage(imageUrl);
        const filename = path.basename(imageUrl);
        const uploadedImage = await uploadImageToSanity(imageBuffer, filename);
        imageMap[imageUrl] = uploadedImage.url;
      } catch (error) {
        console.error(`Failed to process image: ${imageUrl}`, error);
      }
    }

    const updatedContent = await replaceImageUrls(post.content, imageMap);


    const sanityPost = {
      _type: 'blogPostDev',
      title: post.title,
      meta: { slug: { current: post.slug } },
      publishedAt: new Date(post.date).toISOString(),
      body: [
        {
          _type: 'block',
          _key: nanoid(),
          children: [{ _type: 'span', text: updatedContent, _key: nanoid() }],
        },
      ],
    };

    try {
      const result = await client.create(sanityPost);
      console.log(`Successfully migrated post: ${result.title}`);
    } catch (error) {
      console.error(`Failed to migrate post: ${post.title}`, error);
    }
  }
}

migratePosts().then(() => console.log('Migration complete'));