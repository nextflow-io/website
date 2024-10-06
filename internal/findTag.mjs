import sanityClient from '@sanity/client';

export const client = sanityClient({
  projectId: 'o2y1bt2g',
  dataset: 'seqera',
  token: process.env.SANITY_TOKEN,
  useCdn: false,
});

async function findTag(slug) {
  let title = slug.replace(/-/g, ' ');
  if (title === 'pipelines') title = 'seqera pipelines';
  const tag = await client.fetch(`*[_type == "tag" && lower(title) == lower($title)][0]`, { title });
  return tag
}

export default findTag;