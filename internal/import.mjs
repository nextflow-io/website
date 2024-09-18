import sanityClient from '@sanity/client';
import fs from 'fs';
import path from 'path';
import { customAlphabet } from 'nanoid';
import { marked } from 'marked';

const nanoid = customAlphabet('0123456789abcdef', 12);

const client = sanityClient({
  projectId: 'o2y1bt2g',
  dataset: 'seqera',
  token: process.env.SANITY_TOKEN,
  useCdn: false,
});

const postsFile = path.join(process.cwd(), 'export.json');
const contentRoot = path.join(process.cwd(), '../public');

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

function markdownToPortableText(markdown, imageMap) {
  const tokens = marked.lexer(markdown);
  return tokens.map(tokenToPortableText.bind(null, imageMap)).filter(Boolean);
}

function tokenToPortableText(imageMap, token) {
  
  switch (token.type) {
    case 'heading':
      return {
        _type: 'block',
        _key: nanoid(),
        style: `h${token.depth}`,
        children: [{ _type: 'span', text: token.text, _key: nanoid() }],
      };
    case 'paragraph':
      return {
        _type: 'block',
        _key: nanoid(),
        children: token.tokens.map(inlineTokenToPortableText.bind(null, imageMap)),
      };
    case 'image':
      const image = imageMap[src];
      if (!image?._id) {
        console.warn(`Failed to find image for token: ${token.href}`);
        return null;
      }
      return {
        _type: 'image',
        _key: nanoid(),
        asset: {
          _type: 'reference',
          _ref: image._id,
        },
        alt: token.text,
      };
    case 'code':
      return {
        _type: 'code',
        _key: nanoid(),
        code: token.text
      };
    case 'html':
      if (token.text.includes('<img')) {
        const imgTag = token.text.match(/<img.*?>/)[0];
        const srcMatch = imgTag.match(/src=(['"])(.*?)\1/);
        const altMatch = imgTag.match(/alt=(['"])(.*?)\1/);
        const src = srcMatch ? srcMatch[2] : '';
        const alt = altMatch ? altMatch[2] : '';

        const image = imageMap[src];
        if (!image?._id) {
          console.warn(`Failed to find image for token: ${token.text}`);
          return null;
        }
        
        return {
          _type: 'image',
          _key: nanoid(),
          asset: {
            _type: 'reference',
            _ref: image._id,
          },
          alt,
        };
      } else {
        console.warn(`Unsupported HTML token: ${token.text}`);
        return null;
      }
    default:
      console.warn(`Unsupported token type: ${token.type}`, token);
      return null;
  }
}

function inlineTokenToPortableText(imageMap, token) {
  switch (token.type) {
    case 'text':
      let marks = [];
      if (token.bold) marks.push('strong');
      if (token.italic) marks.push('em');
      return { 
        _type: 'span', 
        text: token.text, 
        marks: marks,
        _key: nanoid() 
      };
    case 'link':
      return {
        _type: 'span',
        _key: nanoid(),
        marks: ['link'],
        text: token.text,
        data: { href: token.href },
      };
    case 'image':
      const image = imageMap[token.href];
      if (!image?._id) {
        console.warn(`Failed to find image for token: ${token.href}`);
        return null;
      }
      return {
        _type: 'image',
        _key: nanoid(),
        asset: {
          _type: 'reference',
          _ref: image._id,
        },
        alt: token.text,
      };
    case 'codespan':
      return {
        _type: 'span',
        _key: nanoid(),
        marks: ['code'],
        text: token.text,
      };
    default:
      console.warn(`Unsupported inline token type: ${token.type}`);
      return { _type: 'span', text: token.raw, _key: nanoid() };
  }
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
        imageMap[imagePath] = uploadedImage;
      } catch (error) {
        console.error(`Failed to process image: ${imagePath}`, error);
      }
    }
    console.log('Image map:', imageMap);
    

    const portableTextContent = markdownToPortableText(post.content, imageMap);

    const sanityPost = {
      _type: 'blogPostDev',
      title: post.title,
      meta: { slug: { current: post.slug } },
      publishedAt: new Date(post.date).toISOString(),
      body: portableTextContent,
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