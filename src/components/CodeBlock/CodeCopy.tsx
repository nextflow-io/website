import React, { useEffect, useState, useCallback } from 'react';

export default function CodeCopyButton() {
  const [tooltipVisible, setTooltipVisible] = useState<boolean>(false);
  const [tooltipPosition, setTooltipPosition] = useState<{ x: number, y: number }>({ x: 0, y: 0 });

  const addCopyButtons = useCallback(() => {
    const codeBlocks = document.querySelectorAll('pre > code');
    
    codeBlocks.forEach((codeBlock) => {
      const parent = codeBlock.parentElement;
      if (parent && !parent.querySelector('.copy-button')) {
        const button = document.createElement('button');
        button.className = 'copy-button';
        button.title = 'Copy to clipboard';
        button.innerHTML = `
          <svg width="16" height="16" aria-hidden="true" data-view-component="true">
            <path fill="currentColor" d="M0 6.75C0 5.784.784 5 1.75 5h1.5a.75.75 0 0 1 0 1.5h-1.5a.25.25 0 0 0-.25.25v7.5c0 .138.112.25.25.25h7.5a.25.25 0 0 0 .25-.25v-1.5a.75.75 0 0 1 1.5 0v1.5A1.75 1.75 0 0 1 9.25 16h-7.5A1.75 1.75 0 0 1 0 14.25Z"/>
            <path fill="currentColor" d="M5 1.75C5 .784 5.784 0 6.75 0h7.5C15.216 0 16 .784 16 1.75v7.5A1.75 1.75 0 0 1 14.25 11h-7.5A1.75 1.75 0 0 1 5 9.25Zm1.75-.25a.25.25 0 0 0-.25.25v7.5c0 .138.112.25.25.25h7.5a.25.25 0 0 0 .25-.25v-7.5a.25.25 0 0 0-.25-.25Z"/>
          </svg>
          <span class="copy-tooltip">Copied!</span>
        `;
        
        parent.style.position = 'relative';
        
        button.addEventListener('click', (e) => {
          const code = codeBlock.textContent || '';
          
          navigator.clipboard.writeText(code).then(() => {
            const tooltip = button.querySelector('.copy-tooltip');
            if (tooltip) {
              tooltip.classList.add('visible');
              setTimeout(() => {
                tooltip.classList.remove('visible');
              }, 2000);
            }
          }).catch(err => {
            console.error('Error al copiar: ', err);
          });
        });
        
        parent.appendChild(button);
      }
    });
  }, []);

  const addLinkClassToAnchors = useCallback(() => {
    const paragraphs = document.querySelectorAll('.code-examples p, .blg-summary p');
    
    paragraphs.forEach((paragraph) => {
      const links = paragraph.querySelectorAll('a');
      
      links.forEach((link) => {
        if (!link.classList.contains('link')) {
          link.classList.add('link');
        }
      });
    });
  }, []);

  useEffect(() => {
    const setupPage = () => {
      addCopyButtons();
      addLinkClassToAnchors();
    };

    setupPage();

    const observer = new MutationObserver(setupPage);
    
    observer.observe(document.documentElement, {
      childList: true,
      subtree: true
    });

    window.addEventListener('load', setupPage);

    return () => {
      observer.disconnect();
      window.removeEventListener('load', setupPage);
      document.querySelectorAll('.copy-button').forEach((button) => {
        button.remove();
      });
    };
  }, [addCopyButtons, addLinkClassToAnchors]);

  return null;
} 