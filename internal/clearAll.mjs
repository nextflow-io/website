import sanityClient from '@sanity/client';

export const client = sanityClient({
  projectId: 'o2y1bt2g',
  dataset: 'seqera',
  token: process.env.SANITY_TOKEN,
  useCdn: false,
});

async function deleteAllBlogPostDev() {
  try {
    // 1. Fetch all documents of type 'blogPostDev'
    const query = '*[_type == "blogPostDev"]._id';
    const ids = await client.fetch(query);

    console.log(`Found ${ids.length} blogPostDev documents to delete.`);

    // 2. Delete the documents in batches
    const batchSize = 100; // Adjust based on your needs
    for (let i = 0; i < ids.length; i += batchSize) {
      const batch = ids.slice(i, i + batchSize);
      const transaction = client.transaction();

      batch.forEach(id => {
        transaction.delete(id);
      });

      console.log(`Deleting batch ${i / batchSize + 1}...`);
      await transaction.commit();
      console.log(`Batch ${i / batchSize + 1} deleted.`);
    }

    console.log('All blogPostDev documents have been deleted.');
  } catch (error) {
    console.error('Error deleting documents:', error);
  }
}

deleteAllBlogPostDev();