import React from 'react';
import SideNavigation from './index';

// Ejemplos de uso del componente SideNavigation

export const AmbassadorNav = () => {
  const items = [
    {
      id: 'nextflow-ambassadors',
      title: 'Nextflow ambassadors',
      href: 'nextflow-ambassadors'
    },
    {
      id: 'become-ambassador',
      title: 'Become an ambassador',
      href: 'become-ambassador'
    },
    {
      id: 'our-ambassadors',
      title: 'Our ambassadors',
      href: 'our-ambassadors'
    },
    {
      id: 'fundings',
      title: 'Fundings',
      href: 'fundings'
    }
  ];

  return (
    <SideNavigation 
      items={items} 
      title="Ambassador program"
      activeItem="nextflow-ambassadors" 
    />
  );
};

export const AboutNextflowNav = () => {
  const items = [
    {
      id: 'introduction',
      title: 'Introduction',
      href: 'introduction'
    },
    {
      id: 'key-features',
      title: 'Key features',
      href: 'key-features'
    },
    {
      id: 'community',
      title: 'Community',
      href: 'community'
    },
    {
      id: 'getting-started',
      title: 'Getting started',
      href: 'getting-started'
    },
    {
      id: 'use-cases',
      title: 'Use cases',
      href: 'use-cases'
    }
  ];

  return (
    <SideNavigation 
      items={items} 
      title="About Nextflow"
      activeItem="introduction" 
    />
  );
};

export const ExamplesNav = () => {
  const items = [
    {
      id: 'basic-workflows',
      title: 'Basic workflows',
      href: 'basic-workflows'
    },
    {
      id: 'data-processing',
      title: 'Data processing',
      href: 'data-processing'
    },
    {
      id: 'advanced-patterns',
      title: 'Advanced patterns',
      href: 'advanced-patterns'
    },
    {
      id: 'containers',
      title: 'Containers',
      href: 'containers'
    },
    {
      id: 'nf-core-pipelines',
      title: 'nf-core pipelines',
      href: 'nf-core-pipelines'
    }
  ];

  return (
    <SideNavigation 
      items={items} 
      title="Examples"
      activeItem="basic-workflows" 
    />
  );
}; 