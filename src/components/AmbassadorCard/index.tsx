import React from 'react';
import clsx from 'clsx';
import styles from './styles.module.css';

type AmbassadorCardProps = {
  name: string;
  role: string;
  bio: string;
  avatarUrl: string;
  className?: string;
};

const AmbassadorCard: React.FC<AmbassadorCardProps> = ({
  name,
  role,
  bio,
  avatarUrl,
  className,
}) => {
  return (
    <div className={clsx(styles.card, className)}>
      <div className={styles.avatarContainer}>
        <img
          src={avatarUrl}
          alt={`${name} avatar`}
          className={styles.avatar}
        />
      </div>
      
      <h3 className={styles.name}>{name}</h3>
      <span className={styles.role}>{role}</span>
      
      <p className={styles.bio}>{bio}</p>
    </div>
  );
};

export default AmbassadorCard; 