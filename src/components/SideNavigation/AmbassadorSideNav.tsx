import React from 'react';
import SideNavigation from './index';

const AmbassadorSideNav = () => {
  const items = [
    {
      id: 'introduction',
      title: 'Become an ambassador',
      href: 'introduction'
    },
    {
      id: 'why-become',
      title: 'Why become an ambassador',
      href: 'why-become'
    },
    {
      id: 'what-we-expect',
      title: 'What we expect from you',
      href: 'what-we-expect'
    },
    {
      id: 'become-ambassador',
      title: 'Become an ambassador',
      href: 'become-ambassador'
    },
  ];

  return (
    <SideNavigation 
      items={items} 
      title="Ambassador program"
      activeItem="Become an ambassador" 
      className="font-sans text-sm"
      mode="anchor"
    />
  );
};

export default AmbassadorSideNav; 