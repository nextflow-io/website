import React, { useEffect, useState } from 'react';
import SideNavigation from './index';

const ExamplesSideNav = () => {
  const [activeId, setActiveId] = useState('example1');
  
  // Detectar la página actual en el cliente
  useEffect(() => {
    if (typeof window !== 'undefined') {
      const path = window.location.pathname;
      const filename = path.split('/').pop();
      
      if (filename) {
        // Buscar coincidencia con alguno de los ejemplos
        const match = filename.match(/example(\d+)/);
        if (match && match[1]) {
          setActiveId(`example${match[1]}`);
        }
      }
    }
  }, []);
  
  const items = [
    {
      id: 'example1',
      title: 'Basic pipeline',
      href: 'example1.html'
    },
    {
      id: 'example2',
      title: 'Mixing scripting languages',
      href: 'example2.html'
    },
    {
      id: 'example3',
      title: 'BLAST pipeline',
      href: 'example3.html'
    },
    {
      id: 'example4',
      title: 'RNA-Seq pipeline',
      href: 'example4.html'
    },
    {
      id: 'example5',
      title: 'Machine Learning pipeline',
      href: 'example5.html'
    }
  ];

  return (
    <SideNavigation 
      items={items} 
      title="Examples"
      activeItem={activeId} 
      className="font-sans text-sm"
      mode="page"
    />
  );
};

export default ExamplesSideNav; 