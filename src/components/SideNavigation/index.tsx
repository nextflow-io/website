import React, { useEffect, useRef, useState } from 'react';
import clsx from 'clsx';
import styles from './styles.module.css';

interface NavItem {
  id: string;
  title: string;
  href: string;
}

interface SideNavigationProps {
  items: NavItem[];
  activeItem?: string;
  title?: string;
  className?: string;
  mode?: 'page' | 'anchor'; 
}

interface VisibleSection {
  id: string;
  visiblePercentage: number;
  top: number;
}

function SideNavigation({ items, activeItem, title, className, mode = 'page' }: SideNavigationProps) {
  const [activeId, setActiveId] = useState<string>('');
  const scrollRef = useRef<HTMLDivElement>(null);
  const [isMobile, setIsMobile] = useState<boolean>(false);
  
  useEffect(() => {
    if (typeof window !== 'undefined') {
      const checkIsMobile = () => {
        setIsMobile(window.innerWidth < 1025);
      };
      
      checkIsMobile();
      window.addEventListener('resize', checkIsMobile);
      
      return () => {
        window.removeEventListener('resize', checkIsMobile);
      };
    }
  }, []);
  
  useEffect(() => {
    if (typeof window !== 'undefined') {
      if (mode === 'page') {
        const currentPath = window.location.pathname;
        const currentPage = currentPath.split('/').pop()?.replace('.html', '') || '';
        
        const matchingItem = items.find(item => {
          const itemPath = item.href.replace('.html', '');
          return itemPath === currentPage;
        });
        
        if (matchingItem) {
          setActiveId(matchingItem.id);
        } else if (activeItem) {
          setActiveId(activeItem);
        }
      } else if (mode === 'anchor') {
        const hash = window.location.hash.replace('#', '');
        if (hash) {
          const matchingItem = items.find(item => item.href === hash);
          if (matchingItem) {
            setActiveId(matchingItem.id);
          }
        } else if (activeItem) {
          setActiveId(activeItem);
        }
        
        const handleScroll = () => {
          if (window.isScrolling) return;
          
          const viewportHeight = window.innerHeight;
          const scrollTop = window.scrollY || document.documentElement.scrollTop;
          const sectionIds = items.map(item => item.href);
          
          const visibleSections: VisibleSection[] = sectionIds
            .map(id => {
              const element = document.getElementById(id);
              if (!element) return null;
              
              const rect = element.getBoundingClientRect();
              const elementTop = rect.top + scrollTop;
              const elementBottom = elementTop + element.offsetHeight;
              
              const visibleTop = Math.max(elementTop, scrollTop);
              const visibleBottom = Math.min(elementBottom, scrollTop + viewportHeight);
              const visibleHeight = Math.max(0, visibleBottom - visibleTop);
              
              const percentInView = visibleHeight / element.offsetHeight;
              let bonus = 0;
              
              if (rect.top >= 0 && rect.top < viewportHeight * 0.25) {
                bonus = 0.5; 
              }
              
              return {
                id,
                visiblePercentage: percentInView + bonus,
                top: rect.top
              };
            })
            .filter((section): section is VisibleSection => section !== null);
          
          visibleSections.sort((a, b) => {
            if (Math.abs(a.visiblePercentage - b.visiblePercentage) > 0.2) {
              return b.visiblePercentage - a.visiblePercentage;
            }
            return a.top - b.top;
          });
          
          if (visibleSections.length > 0) {
            const mostVisibleSection = visibleSections[0].id;
            const matchingItem = items.find(item => item.href === mostVisibleSection);
            if (matchingItem && matchingItem.id !== activeId) {
              setActiveId(matchingItem.id);
              
              if (isMobile && scrollRef.current) {
                const activeElement = scrollRef.current.querySelector(`[href="#${matchingItem.href}"]`);
                if (activeElement) {
                  scrollIntoViewIfNeeded(activeElement as HTMLElement, scrollRef.current);
                }
              }
            }
          }
        };
        
        window.addEventListener('scroll', handleScroll, { passive: true });
        
        handleScroll();
        
        return () => {
          window.removeEventListener('scroll', handleScroll);
        };
      }
    }
  }, [items, activeItem, mode, isMobile]);
  
  const scrollIntoViewIfNeeded = (element: HTMLElement, container: HTMLElement) => {
    const containerRect = container.getBoundingClientRect();
    const elementRect = element.getBoundingClientRect();
    
    if (elementRect.left < containerRect.left) {
      container.scrollLeft += elementRect.left - containerRect.left - 20;
    } else if (elementRect.right > containerRect.right) {
      container.scrollLeft += elementRect.right - containerRect.right + 20;
    }
  };
  
  const handleClick = (e: React.MouseEvent<HTMLAnchorElement>, item: NavItem) => {
    if (mode === 'page') {
      if (activeId === item.id) {
        e.preventDefault(); 
      }
    } else if (mode === 'anchor') {
      e.preventDefault(); 
      
      window.isScrolling = true;
      
      const element = document.getElementById(item.href);
      if (element) {
        element.scrollIntoView({ behavior: 'smooth' });
        
        setActiveId(item.id);
        
        if (history.pushState) {
          history.pushState(null, '', `#${item.href}`);
        } else {
          window.location.hash = item.href;
        }
        
        if (isMobile && scrollRef.current) {
          const activeElement = e.currentTarget;
          scrollIntoViewIfNeeded(activeElement, scrollRef.current);
        }
        
        setTimeout(() => {
          window.isScrolling = false;
        }, 1000); 
      }
    }
  };
  
  useEffect(() => {
    if (isMobile && scrollRef.current && activeId) {
      const activeElement = scrollRef.current.querySelector(`a.${styles.activeLink}`);
      if (activeElement) {
        scrollIntoViewIfNeeded(activeElement as HTMLElement, scrollRef.current);
      }
    }
  }, [activeId, isMobile]);
  
  return (
    <nav className={clsx(
      'sticky lg:top-24 lg:self-start lg:w-55',
      'sticky top-[64px] z-20 w-full', 
      className
    )}>
      <div 
        ref={scrollRef}
        className={clsx(
          "lg:border-l border-gray-200",
          "lg:overflow-visible",
          styles.scrollContainer,
          "bg-white pt-2" 
        )}
      >
        {title && (
          <div className={clsx(
            "pb-4 pt-2 lg:py-4 font-['Menlo']",
            styles.titleDisplay
          )}>
            <h3 className="text-[16px] my-0">{title}</h3>
          </div>
        )}
        
        <ul className={clsx(
          "flex lg:flex-col p-0 mt-0",
          styles.mobileNavList
        )}>
          {items.map((item) => (
            <li key={item.id} className="relative group min-w-max">
              <a 
                href={mode === 'page' ? item.href : `#${item.href}`}
                onClick={(e) => handleClick(e, item)}
                className={clsx(
                  "block whitespace-nowrap transition-colors duration-200 mx-1 lg:mx-0 lg:py-4 lg:px-4",
                  activeId === item.id 
                    ? styles.activeLink
                    : " hover:bg-nextflow-light-green hover:text-nextflow-green"
                )}
              >
                {item.title}
              </a>
            </li>
          ))}
        </ul>
        
        <div className={styles.scrollIndicator}></div>
      </div>
    </nav>
  );
}

declare global {
  interface Window {
    isScrolling: boolean;
  }
}

export default SideNavigation; 