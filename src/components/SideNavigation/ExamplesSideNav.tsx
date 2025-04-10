import React, { useEffect, useState } from 'react';
import SideNavigation from './index';

const ExamplesSideNav = () => {
  const [activeId, setActiveId] = useState('basic-pipeline');
  
  useEffect(() => {
    if (typeof window !== 'undefined') {
      const path = window.location.pathname;
      const filename = path.split('/').pop();
      
      if (filename) {
        const match = filename.match(/example(\d+)/);
        if (match && match[1]) {
          setActiveId(`example${match[1]}`);
        }
      }
    }
  }, []);
  
  const items = [
    {
      id: 'basic-pipeline',
      title: 'Basic pipeline',
      href: 'basic-pipeline.html'
    },
    {
      id: 'mixing-scripting-languages',
      title: 'Mixing scripting languages',
      href: 'mixing-scripting-languages.html'
    },
    {
      id: 'blast-pipeline',
      title: 'BLAST pipeline',
      href: 'blast-pipeline.html'
    },
    {
      id: 'rna-seq-pipeline',
      title: 'RNA-Seq pipeline',
      href: 'rna-seq-pipeline.html'
    },
    {
      id: 'machine-learning-pipeline',
      title: 'Machine Learning pipeline',
      href: 'machine-learning-pipeline.html'
    }
  ];

  return (
    <SideNavigation 
      items={items} 
      title="All examples"
      activeItem={activeId} 
      className="font-sans text-sm"
      mode="page"
    />
  );
};

export default ExamplesSideNav; 