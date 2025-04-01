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
        
        // Añadir listener para detectar cambios en el scroll
        const handleScroll = () => {
          // Solo procesar si no estamos en una transición de scroll suave
          if (!window.isScrolling) {
            const viewportHeight = window.innerHeight;
            const sectionIds = items.map(item => item.href);
            
            const visibleSections: VisibleSection[] = sectionIds
              .map(id => {
                const element = document.getElementById(id);
                if (!element) return null;
                
                const rect = element.getBoundingClientRect();
                const visibleHeight = Math.min(rect.bottom, viewportHeight) - Math.max(rect.top, 0);
                const visiblePercentage = visibleHeight > 0 ? visibleHeight / element.offsetHeight : 0;
                
                return {
                  id,
                  visiblePercentage,
                  top: rect.top
                };
              })
              .filter((section): section is VisibleSection => section !== null)
              .sort((a, b) => {
                if (a.visiblePercentage !== b.visiblePercentage) {
                  return b.visiblePercentage - a.visiblePercentage;
                }
                return a.top - b.top;
              });
            
            if (visibleSections.length > 0 && visibleSections[0].visiblePercentage > 0) {
              const mostVisibleSection = visibleSections[0].id;
              const matchingItem = items.find(item => item.href === mostVisibleSection);
              if (matchingItem && matchingItem.id !== activeId) {
                setActiveId(matchingItem.id);
              }
            }
          }
        };
        
        window.addEventListener('scroll', handleScroll);
        handleScroll(); // Ejecutar inmediatamente
        
        return () => {
          window.removeEventListener('scroll', handleScroll);
        };
      }
    }
  }, [items, activeItem, mode]);
  
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
        
        setTimeout(() => {
          window.isScrolling = false;
        }, 1000); 
      }
    }
  };
  
  return (
    <nav className={clsx(
      'lg:w-64 lg:sticky lg:top-24 lg:self-start',
      'w-full sticky top-0 z-10',
      className
    )}>
      <div 
        ref={scrollRef}
        className={clsx(
          "lg:border-l border-gray-200",
          "lg:overflow-visible",
          styles.scrollContainer
        )}
      >
        {title && (
          <div className={clsx(
            "py-4 lg:py-0 font-['Menlo']",
            styles.titleDisplay
          )}>
            <h3 className="text-[16px]">{title}</h3>
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
                    : "text-gray-900 hover:bg-nextflow-light-green hover:text-nextflow-green"
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