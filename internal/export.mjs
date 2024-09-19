import fs from 'fs';
import path from 'path';
import matter from 'gray-matter';
import * as cheerio from 'cheerio';

const postsDirectory = path.join(process.cwd(), '../src/content/blog');
const outputFile = path.join(process.cwd(), 'export.json');

function extractImagePaths(content, postPath) {
  const $ = cheerio.load(content);
  const images = [];
  $('img').each((i, elem) => {
    const src = $(elem).attr('src');
    if (src) {
      images.push(src);
    }
  });
  return images;
}

function getPostsRecursively(dir) {
  let posts = [];
  const items = fs.readdirSync(dir, { withFileTypes: true });

  for (const item of items) {
    const fullPath = path.join(dir, item.name);
    
    if (item.isDirectory()) {
      posts = posts.concat(getPostsRecursively(fullPath));
    } else if (item.isFile() && item.name.endsWith('.md')) {
      const fileContents = fs.readFileSync(fullPath, 'utf8');
      const { data, content } = matter(fileContents);
      const images = extractImagePaths(content, fullPath);
      
      posts.push({
        slug: path.relative(postsDirectory, fullPath).replace('.md', ''),
        title: data.title,
        date: data.date,
        content: content,
        images: images,
      });
    }
  }

  return posts;
}

const posts = getPostsRecursively(postsDirectory);

fs.writeFileSync(outputFile, JSON.stringify(posts, null, 2));
console.log(`Exported ${posts.length} posts to ${outputFile}`);