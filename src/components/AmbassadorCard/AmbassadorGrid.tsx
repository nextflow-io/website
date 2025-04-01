import React from 'react';
import clsx from 'clsx';
import AmbassadorCard from './index';
import styles from './styles.module.css';

type Ambassador = {
  id: string;
  name: string;
  img: string;
  country: string;
  github: string;
  linkedin?: string;
  twitter?: string;
  mastodon?: string;
  bluesky?: string;
  bio?: React.ReactNode;
  title?: string;
};

type AmbassadorGridProps = {
  ambassadors: Ambassador[];
  className?: string;
};

const AmbassadorGrid: React.FC<AmbassadorGridProps> = ({ 
  ambassadors,
  className 
}) => {
  return (
    <div className={clsx(styles.grid, className)}>
      {ambassadors.map((ambassador) => (
        <AmbassadorCard
          key={ambassador.id}
          name={ambassador.name}
          img={ambassador.img}
          country={ambassador.country}
          github={ambassador.github}
          linkedin={ambassador.linkedin}
          twitter={ambassador.twitter}
          mastodon={ambassador.mastodon}
          bluesky={ambassador.bluesky}
          title={ambassador.title}
        >
          {ambassador.bio}
        </AmbassadorCard>
      ))}
    </div>
  );
};

export default AmbassadorGrid; 