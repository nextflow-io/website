import sanityClient from '@sanity/client';

export const client = sanityClient({
  projectId: 'o2y1bt2g',
  dataset: 'seqera',
  token: process.env.SANITY_TOKEN,
  useCdn: false,
});

async function findPerson(name) {
  const person = await client.fetch(`*[_type == "person" && name == $name][0]`, { name });
  if (!person) {
    console.log(`⭕ No person found with the name "${name}".`);
    return;
  } else {
    console.log(`Person found`, person.name);
    return person;
  }
}

export default findPerson;