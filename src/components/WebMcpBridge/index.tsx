import { useEffect } from "react";

type WebMcpTool = {
  name: string;
  description: string;
  inputSchema: Record<string, unknown>;
  execute: (input: Record<string, unknown>) => Promise<unknown>;
};

const tools: WebMcpTool[] = [
  {
    name: "get_documentation_links",
    description: "Return canonical Nextflow documentation and training URLs.",
    inputSchema: {
      type: "object",
      properties: {},
    },
    execute: async () => ({
      documentation: "https://docs.seqera.io/nextflow/",
      training: "https://training.nextflow.io/latest/",
      community: "https://community.seqera.io/tag/nextflow",
      examples: "https://nextflow.io/examples.html",
      api_catalog: "https://nextflow.io/.well-known/api-catalog",
    }),
  },
  {
    name: "get_examples",
    description: "List Nextflow pipeline example pages available on nextflow.io.",
    inputSchema: {
      type: "object",
      properties: {},
    },
    execute: async () => [
      { title: "Basic pipeline", url: "https://nextflow.io/basic-pipeline.html" },
      { title: "Mixing scripting languages", url: "https://nextflow.io/mixing-scripting-languages.html" },
      { title: "BLAST pipeline", url: "https://nextflow.io/blast-pipeline.html" },
      { title: "RNA-seq pipeline", url: "https://nextflow.io/rna-seq-pipeline.html" },
      { title: "Machine learning pipeline", url: "https://nextflow.io/machine-learning-pipeline.html" },
    ],
  },
  {
    name: "search_site",
    description: "Search nextflow.io pages by keyword against known site resources.",
    inputSchema: {
      type: "object",
      properties: {
        query: { type: "string", description: "Search keyword" },
      },
      required: ["query"],
    },
    execute: async ({ query }) => {
      const q = String(query ?? "").toLowerCase();
      const pages = [
        { title: "Home", url: "https://nextflow.io/", keywords: ["nextflow", "workflow", "pipeline", "reproducible"] },
        { title: "Examples", url: "https://nextflow.io/examples.html", keywords: ["example", "tutorial", "pipeline"] },
        { title: "About", url: "https://nextflow.io/about-us.html", keywords: ["about", "seqera", "team"] },
        { title: "Basic pipeline", url: "https://nextflow.io/basic-pipeline.html", keywords: ["basic", "hello", "world"] },
        { title: "Documentation", url: "https://docs.seqera.io/nextflow/", keywords: ["docs", "documentation", "reference", "manual"] },
        { title: "Training", url: "https://training.nextflow.io/latest/", keywords: ["training", "course", "learn"] },
      ];
      return pages.filter(
        (page) =>
          page.title.toLowerCase().includes(q) ||
          page.keywords.some((keyword) => keyword.includes(q) || q.includes(keyword)),
      );
    },
  },
];

export default function WebMcpBridge() {
  useEffect(() => {
    const modelContext = (navigator as Navigator & { modelContext?: Record<string, unknown> }).modelContext;
    if (!modelContext) {
      return;
    }

    const abortController = new AbortController();

    if (typeof modelContext.provideContext === "function") {
      (modelContext.provideContext as (ctx: { tools: WebMcpTool[] }) => void)({ tools });
    } else if (typeof modelContext.registerTool === "function") {
      for (const tool of tools) {
        (modelContext.registerTool as (tool: WebMcpTool, opts: { signal: AbortSignal }) => void)(tool, {
          signal: abortController.signal,
        });
      }
    }

    return () => {
      abortController.abort();
    };
  }, []);

  return null;
}
