import React from 'react';
import AmbassadorCard from './index';
import styles from './styles.module.css';

type Ambassador = {
  id: string;
  name: string;
  role: string;
  bio: string;
  avatarUrl: string;
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
    <div className={styles.grid}>
      {ambassadors.map((ambassador) => (
        <AmbassadorCard
          key={ambassador.id}
          name={ambassador.name}
          role={ambassador.role}
          bio={ambassador.bio}
          avatarUrl={ambassador.avatarUrl}
        />
      ))}
    </div>
  );
};

export default AmbassadorGrid; 