import fs from 'fs';
import path from 'path';
import sanityClient from '@sanity/client';

const outputFile = path.join(process.cwd(), './links2.csv');

export const client = sanityClient({
  projectId: 'o2y1bt2g',
  dataset: 'seqera',
  token: process.env.SANITY_TOKEN,
  useCdn: false,
});

async function fetchNewPosts() {
  return await client.fetch(`*[_type == "blogPostDev2"]`);
}

async function getLinks() {
  console.log('🟢🟢🟢 Export');
  const newPosts = await fetchNewPosts();

  let csvContent = 'title,oldURL,devURL,prodURL,cmsURL\n';

  for (const post of newPosts) {

    let id = post?._id || '';
    if (id.split('.')[1]) id = id.split('.')[1];

    const slug = post?.meta?.slug?.current || '';
    const title = post?.title || '';

    const oldURL = post ? `https://seqera.io/blog/${slug}` : '';
    const devURL = post ? `https://seqera.io/preview?type=blogPostDev&id=${id}` : '';
    const prodURL = ''
    const cmsURL = 'https://seqera-cms.netlify.app/seqera/structure/blogPostDev2;' + id;

    const escapedTitle = title.includes(',') ? `"${title}"` : title;
    csvContent += `${escapedTitle},${oldURL},${devURL},${prodURL},${cmsURL}\n`;
  }

  fs.writeFileSync(outputFile, csvContent, 'utf8');
  console.log(`CSV file has been written to ${outputFile}`);
}

getLinks();