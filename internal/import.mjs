import sanityClient from '@sanity/client';
import fs from 'fs';
import path from 'path';
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
const contentRoot = path.join(process.cwd(), '../public'); // Adjust this to your content root

async function readPosts() {
  const data = await fs.promises.readFile(postsFile, 'utf8');
  return JSON.parse(data);
}

async function readImageFromFileSystem(imagePath) {
  const fullPath = path.join(contentRoot, imagePath);
  return fs.promises.readFile(fullPath);
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
  const posts = await readPosts();
  const p = [posts[4]];

  for (const post of p) {
    const imageMap = {};
    for (const imagePath of post.images) {
      try {
        const imageBuffer = await readImageFromFileSystem(imagePath);
        const filename = path.basename(imagePath);
        const uploadedImage = await uploadImageToSanity(imageBuffer, filename);
        imageMap[imagePath] = uploadedImage.url;
      } catch (error) {
        console.error(`Failed to process image: ${imagePath}`, error);
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