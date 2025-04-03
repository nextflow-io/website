import React from 'react';
import SideNavigation from './index';

const AboutSideNav = () => {
  const items = [
    {
      id: 'about-nextflow',
      title: 'About Nextflow',
      href: 'about-nextflow'
    },
    {
      id: 'core-team',
      title: 'Core team',
      href: 'core-team'
    },
    {
      id: 'project-funding',
      title: 'Project Funding',
      href: 'project-funding'
    },
    {
      id: 'key-contributors',
      title: 'Key contributors',
      href: 'key-contributors'
    }
  ];

  return (
    <SideNavigation 
      items={items} 
      title="About"
      activeItem="about-nextflow" 
      className="font-sans text-sm"
      mode="anchor"
    />
  );
};

export default AboutSideNav; 