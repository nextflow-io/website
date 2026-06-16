const MARKDOWN_PATHS = {
  "/": "/index.md",
  "/index.html": "/index.md",
};

export default async function handler(request, context) {
  const accept = request.headers.get("Accept") ?? "";
  if (!accept.includes("text/markdown")) {
    return context.next();
  }

  const url = new URL(request.url);
  const markdownPath = MARKDOWN_PATHS[url.pathname];
  if (!markdownPath) {
    return context.next();
  }

  const markdownUrl = new URL(markdownPath, url.origin);
  const markdownResponse = await context.rewrite(markdownUrl.toString());
  if (!markdownResponse.ok) {
    return context.next();
  }

  const body = await markdownResponse.text();
  const tokenEstimate = Math.ceil(body.length / 4);

  return new Response(body, {
    status: 200,
    headers: {
      "Content-Type": "text/markdown; charset=utf-8",
      "x-markdown-tokens": String(tokenEstimate),
      Vary: "Accept",
      "Cache-Control": "public, max-age=300, must-revalidate",
    },
  });
}
