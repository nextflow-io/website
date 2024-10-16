import fs from 'fs';
import path from 'path';
import sanityClient from '@sanity/client';

const outputFile = path.join(process.cwd(), 'backup3.json');

export const client = sanityClient({
  projectId: 'o2y1bt2g',
  dataset: 'seqera',
  token: process.env.SANITY_TOKEN,
  useCdn: false,
});


async function fetchBlogPosts() {
  return await client.fetch(`*[_type == "blogPost"]`);
}

async function doBackup() {
  const posts = await fetchBlogPosts();
  fs.writeFileSync(outputFile, JSON.stringify(posts, null, 2));
}

doBackup();