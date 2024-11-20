import fs from 'fs';
import path from 'path';
import matter from 'gray-matter';
import * as cheerio from 'cheerio';

const postsDirectory = path.join(process.cwd(), '../src/content/blog');
const outputFile = path.join(process.cwd(), 'export.json');

function extractImagePaths(content, postPath) {
  const images = [];

  // Extract HTML images
  const $ = cheerio.load(content);
  $('img').each((i, elem) => {
    const src = $(elem).attr('src');
    if (src) {
      images.push(src);
    }
  });

  // Extract Markdown images
  const markdownImageRegex = /!\[.*?\]\((.*?)\)/g;
  let match;
  while ((match = markdownImageRegex.exec(content)) !== null) {
    images.push(match[1]);
  }

  // Remove duplicates
  return [...new Set(images)];
}

function sanitizeMarkdown(content) {
  const $ = cheerio.load(`<div id="root">${content}</div>`);

  $('p').each((i, elem) => {
    const $elem = $(elem);
    $elem.replaceWith(`\n\n${$elem.html().trim()}\n\n`);
  });

  $('s, del, strike').each((i, elem) => {
    const $elem = $(elem);
    $elem.replaceWith(`~~${$elem.html()}~~`);
  });

  $('sup').each((i, elem) => {
    const $elem = $(elem);
    $elem.replaceWith(`^${$elem.html()}^`);
  });

  $('sub').each((i, elem) => {
    const $elem = $(elem);
    $elem.replaceWith(`~${$elem.html()}~`);
  });

  $('a').each((i, elem) => {
    const $elem = $(elem);
    const href = $elem.attr('href');
    const text = $elem.text().replace(/\n/g, ' ');
    $elem.replaceWith(`[${text}](${href})`);
  });

  $('blockquote').each((i, elem) => {
    const $elem = $(elem);
    const text = $elem.html().trim().replace(/\n/g, '\n> ');
    $elem.replaceWith(`\n\n> ${text}\n\n`);
  });

  $('em, i').each((i, elem) => {
    const $elem = $(elem);
    $elem.replaceWith(`*${$elem.html().replace(/\n/g, ' ')}*`);
  });

  $('strong, b').each((i, elem) => {
    const $elem = $(elem);
    $elem.replaceWith(`**${$elem.html().replace(/\n/g, ' ')}**`);
  });

  $('code').each((i, elem) => {
    const $elem = $(elem);
    if ($elem.parent().is('pre')) {
      // This is a code block, leave it as is
      return;
    }
    $elem.replaceWith(`\`${$elem.html()}\``);
  });

  $('hr').each((i, elem) => {
    $(elem).replaceWith('\n\n---\n\n');
  });

  $('ul, ol').each((i, elem) => {
    const $elem = $(elem);
    const listItems = $elem.children('li').map((i, li) => {
      const prefix = $elem.is('ul') ? '- ' : `${i + 1}. `;
      return prefix + $(li).html().trim();
    }).get().join('\n');
    $elem.replaceWith(`\n\n${listItems}\n\n`);
  });

  // Remove any remaining HTML tags
  // $('*').each((i, elem) => {
  //   const $elem = $(elem);
  //   $elem.replaceWith($elem.html());
  // });

  let markdown = $('#root').html().trim();
  markdown = markdown.replace(/\n{3,}/g, '\n\n');
  return markdown;
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
      const convertedContent = sanitizeMarkdown(content);
      const images = extractImagePaths(convertedContent, fullPath);

      posts.push({
        slug: path.relative(postsDirectory, fullPath).replace('.md', ''),
        title: data.title,
        date: data.date,
        content: convertedContent,
        images: images,
        author: data.author,
        tags: data.tags,
      });
    }
  }

  return posts;
}

const posts = getPostsRecursively(postsDirectory);

fs.writeFileSync(outputFile, JSON.stringify(posts, null, 2));
console.log(`Exported ${posts.length} posts to ${outputFile}`);