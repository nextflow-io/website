import React, { useEffect, useState } from 'react';
import SideNavigation from './index';

const ExamplesSideNav = () => {
  const [activeId, setActiveId] = useState('basic-pipeline');
  
  useEffect(() => {
    if (typeof window !== 'undefined') {
      const path = window.location.pathname;

      // Extract pipeline name from /examples/pipeline-name/ format
      const match = path.match(/\/examples\/([^\/]+)\/?$/);
      if (match && match[1]) {
        setActiveId(match[1]);
      }
    }
  }, []);
  
  const items = [
    {
      id: 'basic-pipeline',
      title: 'Basic pipeline',
      href: '/examples/basic-pipeline'
    },
    {
      id: 'mixing-scripting-languages',
      title: 'Mixing scripting languages',
      href: '/examples/mixing-scripting-languages'
    },
    {
      id: 'blast-pipeline',
      title: 'BLAST pipeline',
      href: '/examples/blast-pipeline'
    },
    {
      id: 'rna-seq-pipeline',
      title: 'RNA-Seq pipeline',
      href: '/examples/rna-seq-pipeline'
    },
    {
      id: 'machine-learning-pipeline',
      title: 'Machine Learning pipeline',
      href: '/examples/machine-learning-pipeline'
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