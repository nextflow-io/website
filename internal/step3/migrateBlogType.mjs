import sanityClient from '@sanity/client';
import { customAlphabet } from 'nanoid';

const nanoid = customAlphabet('0123456789abcdef', 12);

export const client = sanityClient({
  projectId: 'o2y1bt2g',
  dataset: 'seqera',
  token: process.env.SANITY_TOKEN,
  useCdn: false,
});

async function fetchBlogPostsDev() {
  return await client.fetch(`*[_type == "blogPostDev"]`);
}

async function fetchBlogPosts() {
  return await client.fetch(`*[_type == "blogPost"]`);
}

async function migrateBlogType() {
  console.log('🟢🟢🟢 Migrating');
  const devPosts = await fetchBlogPostsDev();
  const posts = await fetchBlogPosts();

  for (const post of devPosts) {
    console.log('🔵 >>         ', post.meta.slug.current);
    const existing = posts.find(p => p.meta.slug.current === post.meta.slug.current);
    if (!!existing) {
      console.log('🟡 exists   >> ', existing.meta.slug.current);
      console.log('🟡 skipping >> ', existing.title);
      continue;
    }
    const newPost = {
      ...post,
      _type: 'blogPost',
      _id: nanoid(),
      _rev: undefined,
    }
    const p = await client.createIfNotExists(newPost);
    console.log('🟢 created >> ', p.title);
  }
  console.log('🟢🟢🟢 Done');
}

migrateBlogType();