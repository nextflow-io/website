import sanityClient from '@sanity/client';
import fs from 'fs';
import path from 'path';
import axios from 'axios';
import cheerio from 'cheerio';

const client = sanityClient({
  projectId: 'your-project-id',
  dataset: 'your-dataset',
  token: 'your-write-token',
  useCdn: false,
});

const postsFile = path.join(process.cwd(), 'blog-posts.json');
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
  for (const post of posts) {
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
      _type: 'post',
      title: post.title,
      slug: { current: post.slug },
      publishedAt: new Date(post.date).toISOString(),
      body: [
        {
          _type: 'block',
          children: [{ _type: 'span', text: updatedContent }],
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