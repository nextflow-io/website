import rss from "@astrojs/rss";
import { getCollection } from "astro:content";

export async function GET(context) {
  const blog = (await getCollection("blog")).sort(
    (a, b) => new Date(b.data.date).valueOf() - new Date(a.data.date).valueOf(),
  );
  return rss({
    title: "Nextflow Blog",
    description: "Blogging about Nextflow, computational pipelines and parallel programming",
    site: context.site,
    items: blog.map((post) => ({
      title: post.data.title,
      pubDate: post.data.date,
      description: post.data.description,
      link: `/blog/${post.slug}.html`,
    })),
    customData: `<language>en-gb</language>`,
  });
}
