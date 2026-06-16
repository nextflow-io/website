// CloudFront Function: serve markdown when Accept: text/markdown is present.
// Associate with viewer-request on the nextflow.io CloudFront distribution.
function handler(event) {
  var request = event.request;
  var headers = request.headers;
  var accept = headers.accept ? headers.accept.value : "";
  var uri = request.uri;

  if (accept.indexOf("text/markdown") === -1) {
    return request;
  }

  if (uri === "/" || uri === "/index.html") {
    request.uri = "/index.md";
  }

  return request;
}
