import React from 'react';
import SideNavigation from './index';

const ExamplesSideNav = () => {
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
      activeItem="basic-pipeline"
      className="font-sans text-sm"
      mode="page"
    />
  );
};

export default ExamplesSideNav;
