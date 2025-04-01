import React from 'react';
import SideNavigation from './index';

const AmbassadorSideNav = () => {
  const items = [
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
    {
      id: 'our-ambassadors',
      title: 'Our ambassadors',
      href: 'our-ambassadors'
    }
  ];

  return (
    <SideNavigation 
      items={items} 
      title="Ambassador program"
      activeItem="why-become" 
      className="font-sans text-sm"
      mode="anchor"
    />
  );
};

export default AmbassadorSideNav; 