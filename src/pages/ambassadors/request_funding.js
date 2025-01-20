export async function get() {
  return new Response(null, {
    status: 302,
    headers: {
      Location: 'https://script.google.com/macros/s/AKfycbwCXBhHS6mE7Zl-PIoRvMplloyaflhLPVLy82nUyYvC-i9x_iAQ2Cdjtr99L0VsPTE/exec',
    },
  });
}

