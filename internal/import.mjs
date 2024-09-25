import sanityClient from '@sanity/client';
import fs from 'fs';
import path from 'path';
import { customAlphabet } from 'nanoid';
import { marked } from 'marked';
import findPerson from './findPerson.mjs';

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

function sanitizeText(text, removeLineBreaks = false) {
  // Replace all instances of &#39; with '
  const t = text.replace(/&#39;/g, "'");

  if (removeLineBreaks) return t.replace(/\n/g, ' ');

  return t
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
      const children = [];
      const markDefs = [];

      token.tokens.forEach(t => {
        if (t.type === 'link') {
          const linkKey = nanoid();
          children.push({
            _type: 'span',
            _key: nanoid(),
            marks: [linkKey],
            text: sanitizeText(t.text)
          });
          markDefs.push({
            _key: linkKey,
            _type: 'link',
            href: t.href
          });
        } else {
          children.push(inlineTokenToPortableText(imageMap, t));
        }
      });

      return {
        _type: 'block',
        _key: nanoid(),
        style: 'normal',
        children: children,
        markDefs: markDefs
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
        listItem: 'bullet',
        style: 'normal',
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
        text: sanitizeText(token.text, true),
        marks: marks,
        _key: nanoid()
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
        text: sanitizeText(token.text, true),
        marks: [],
      };
    default:
      console.warn(`Unsupported inline token type: ${token.type}`, token);
      return { _type: 'span', text: token.raw, _key: nanoid() };
  }
}

async function migratePosts() {
  const posts = await readPosts();
  const firstTen = posts.slice(0, 10);
  const selected = [
    '2016/deploy-in-the-cloud-at-snap-of-a-finger',
    '2017/caw-and-singularity',
  ]
  const selectedPosts = posts.filter(post => selected.includes(post.slug));

  console.log('');
  console.log('');
  console.log('');
  console.log('');
  console.log('');
  console.log('');
  console.log('');
  console.log('');
  console.log('');
  console.log('🟢🟢🟢 Migrating posts...');
  console.log('');


  for (const post of selectedPosts) {

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

    const person = await findPerson(post.author);
    if (!person) return false;

    const portableTextContent = markdownToPortableText(post.content, imageMap);

    const newSlug = post.slug.split('/').pop();

    let dateStr = post.date.split('T')[0];
    dateStr = `${dateStr} 8:00`;
    console.log(dateStr);


    const sanityPost = {
      _type: 'blogPostDev',
      title: post.title,
      meta: { slug: { current: newSlug } },
      publishedAt: new Date(dateStr).toISOString(),
      body: portableTextContent,
      author: { _type: 'reference', _ref: person._id },
    };

    try {
      const result = await client.create(sanityPost);
      console.log(`✅ Successfully migrated post: ${result.title}`);
    } catch (error) {
      console.error(`Failed to migrate post: ${post.title}`, error);
    }
  }
}

migratePosts().then((isSuccess) => {
  if (isSuccess) console.log('✅ Migration complete')
  else console.log('❌ Migration failed')
});