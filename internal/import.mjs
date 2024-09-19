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

function sanitizeText(text) {
  // Replace all instances of &#39; with '
  return text.replace(/&#39;/g, "'");
}

function tokenToPortableText(imageMap, token) {
  
  switch (token.type) {
    case 'heading':
      return {
        _type: 'block',
        _key: nanoid(),
        style: `h${token.depth}`,
        children: [{ _type: 'span', text: sanitizeText(token.text), _key: nanoid() }],
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

      } else if (token.text.startsWith('<script')) {

        const idMatch = token.text.match(/id=(['"])(.*?)\1/);
        const srcMatch = token.text.match(/src=(['"])(.*?)\1/);
        const id = idMatch ? idMatch[2] : '';
        const src = srcMatch ? srcMatch[2] : '';
        if (!src) {
          console.warn(`Failed to find src for script: ${token.text}`);
          return null;
        }
        return { _type: 'script', _key: nanoid(), id, src };

      } else {
        console.warn(`Unsupported HTML token: ${token.text}`);
        return null;
      }
    case 'list':
      return {
        _type: 'block',
        _key: nanoid(),
        style: token.ordered ? 'number' : 'bullet',
        children: token.items.flatMap(item => 
          item.tokens.map(inlineTokenToPortableText.bind(null, imageMap))
        ),
      };
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
        text: sanitizeText(token.text), 
        marks: marks,
        _key: nanoid() 
      };
    case 'link':
      return {
        _type: 'span',
        _key: nanoid(),
        marks: ['link'],
        text: sanitizeText(token.text),
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
    case 'list_item':
      return {
        _type: 'span',
        _key: nanoid(),
        text: sanitizeText(token.text),
      };
    default:
      console.warn(`Unsupported inline token type: ${token.type}`);
      return { _type: 'span', text: token.raw, _key: nanoid() };
  }
}

async function migratePosts() {
  const posts = await readPosts();
  const firstTen = posts.slice(0, 10);
  const selectedPost = posts.find(p => p.slug === '2016/deploy-in-the-cloud-at-snap-of-a-finger');

  console.log('');
  console.log('');
  console.log('');
  console.log('');
  console.log('');
  console.log('');
  console.log('🪣 Migrating posts...');
  console.log('');
  

  for (const post of [selectedPost]) {
    
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
    
    const portableTextContent = markdownToPortableText(post.content, imageMap);

    const newSlug = post.slug.split('/').pop();

    const sanityPost = {
      _type: 'blogPostDev',
      title: post.title,
      meta: { slug: { current: newSlug } },
      publishedAt: new Date(post.date).toISOString(),
      body: portableTextContent,
    };

    try {
      const result = await client.create(sanityPost);
      console.log(`✅ Successfully migrated post: ${result.title}`);
    } catch (error) {
      console.error(`Failed to migrate post: ${post.title}`, error);
    }
  }
}

migratePosts().then(() => console.log('Migration complete'));