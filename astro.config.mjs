import { defineConfig } from 'astro/config';

import remarkDescription from 'astro-remark-description';

// https://astro.build/config
export default defineConfig({
    site: 'https://nextflow.io/',
    markdown: {
        remarkPlugins: [
            [
                remarkDescription,
                {
                    name: 'excerpt',
                    node: (node, i, parent) => {
                        // check if parent has a child that is an html comment with the text 'end-archive-description'
                        if (
                            parent?.children?.some(
                                (child) =>
                                    (child.type === 'html' && child.value === '<!-- end-archive-description -->') ||
                                    (child.type === 'mdxFlowExpression' && child?.value === '/* end-archive-description */'),
                            )
                        ) {
                            const sibling = parent?.children[i + 1];

                            return (
                                (sibling?.type === 'html' && sibling?.value === '<!-- end-archive-description -->') ||
                                (sibling?.type === 'mdxFlowExpression' && sibling?.value === '/* end-archive-description */')
                            );
                        } else {
                            // return the first paragraph otherwise

                            // get the index of the first paragraph
                            const firstParagraphIndex = parent?.children.findIndex(
                                (child) => child.type === 'paragraph',
                            );
                            // if the node is the first paragraph, return true
                            return i === firstParagraphIndex;
                        }
                    },
                },
            ]
        ]
    }
});
