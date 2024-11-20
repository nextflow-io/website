import fs from 'fs';
import path from 'path';

const postsFile = path.join(process.cwd(), '../export.json');
const outputFile = path.join(process.cwd(), 'redirects.txt');

async function readPosts() {
  const data = await fs.promises.readFile(postsFile, 'utf8');
  return JSON.parse(data);
}

async function getRedirects() {
  const posts = await readPosts();
  let txt = '';
  for (const post of posts) {
    const newSlug = post.slug.split('/').pop();
    txt += `/blog/${post.slug}.html  https://seqera.io/blog/${newSlug}  301`;
    txt += "\n";
  }
  fs.writeFileSync(outputFile, txt);
}

getRedirects();