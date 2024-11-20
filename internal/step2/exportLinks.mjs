import fs from 'fs';
import path from 'path';
import sanityClient from '@sanity/client';

const titlesFile = path.join(process.cwd(), 'posts.csv');
const postsFile = path.join(process.cwd(), '../export.json');
const outputFile = path.join(process.cwd(), './links.csv');

async function readPosts() {
  const data = await fs.promises.readFile(postsFile, 'utf8');
  return JSON.parse(data);
}

export const client = sanityClient({
  projectId: 'o2y1bt2g',
  dataset: 'seqera',
  token: process.env.SANITY_TOKEN,
  useCdn: false,
});

async function fetchNewPosts() {
  return await client.fetch(`*[_type == "blogPostDev"]`);
}

async function getLinks() {
  console.log('🟢🟢🟢 Export');
  const fileContents = fs.readFileSync(titlesFile, 'utf8');
  const titles = fileContents.split('\n');

  const oldPosts = await readPosts();
  const newPosts = await fetchNewPosts();

  let csvContent = 'title,oldURL,devURL,prodURL,cmsURL\n';

  for (const title of titles) {
    const oldPost = oldPosts.find(p => p.title === title);
    const newPost = newPosts.find(p => p.title === title);
    if (!oldPost) console.log('⭕ old: ', title);
    if (!newPost) console.log('⭕ new: ', title);

    let id = newPost?._id || '';
    if (id.split('.')[1]) id = id.split('.')[1];

    let newSlug = newPost?.meta?.slug?.current || '';

    const oldURL = oldPost ? `https://nextflow.io/blog/${oldPost.slug}.html` : '';
    const devURL = newPost ? `https://seqera.io/preview?type=blogPostDev&id=${id}` : '';
    const prodURL = newPost ? `https://seqera.io/blog/${newSlug}` : '';
    const cmsURL = 'https://seqera-cms.netlify.app/seqera/structure/blogPostDev;' + id;

    const escapedTitle = title.includes(',') ? `"${title}"` : title;
    csvContent += `${escapedTitle},${oldURL},${devURL},${prodURL},${cmsURL}\n`;
  }

  fs.writeFileSync(outputFile, csvContent, 'utf8');
  console.log(`CSV file has been written to ${outputFile}`);
}

getLinks();