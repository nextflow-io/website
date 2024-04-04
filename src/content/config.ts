// 1. Import utilities from `astro:content`
import { defineCollection, z } from "astro:content";

const blogCollection = defineCollection({
  type: "content",
  schema: z.object({
    title: z.string(),
    description: z.string().optional(),
    author: z.string(),
    image: z
      .string()
      .regex(/^\/img\//)
      .optional(),
    author2: z.string().optional(),
    icon: z.string(),
    icon2: z.string().optional(),
    community_post: z.boolean().optional(),
    ambassador_post: z.boolean().optional(),
    tags: z.string(),
    // tags: z.array(z.string()).optional(),
    date: z
      .date()
      .transform((d) => d.getDate() + " " + d.toLocaleString("default", { month: "long" }) + " " + d.getFullYear()),
  }),
});

const podcastCollection = defineCollection({
  type: "content",
  schema: z.object({
    title: z.string(),
    episode: z.number(),
    description: z.string(),
    subtype: z.string(),
    author: z.string(),
    youtubeid: z.string().length(11),
    image: z.string().regex(/^\/img\//),
    icon: z.string(),
    tags: z.string(),
    // tags: z.array(z.string()).optional(),
    date: z
      .date()
      .transform((d) => d.getDate() + " " + d.toLocaleString("default", { month: "long" }) + " " + d.getFullYear()),
  }),
});

export const collections = {
  blog: blogCollection,
  podcast: podcastCollection,
};
