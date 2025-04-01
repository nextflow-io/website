import React from 'react';
import clsx from 'clsx';
import styles from './styles.module.css';

export type AmbassadorCardProps = {
  name: string;
  img: string;
  country: string;
  github: string;
  linkedin?: string;
  twitter?: string;
  mastodon?: string;
  bluesky?: string;
  className?: string;
  children?: React.ReactNode;
  title?: string;
};

const AmbassadorCard: React.FC<AmbassadorCardProps> = ({
  name,
  img,
  country,
  github,
  linkedin,
  twitter,
  mastodon,
  bluesky,
  className,
  children,
  title = 'Nextflow Ambassador',
}) => {
  return (
    <div className={clsx(styles.card, className)}>
      <div className={styles.cardImage}>
        <img src={`img/${img}`} className={styles.avatar} alt={`${name}`} />
        
        {country && (
          <div className={styles.flag}>
            <img
              src={`https://flagicons.lipis.dev/flags/4x3/${country}.svg`}
              width="30"
              height="30"
              alt={`${country} flag`}
            />
          </div>
        )}
      </div>
      
      <div className={styles.cardContent}>
        <h2 className={styles.name}>{name}</h2>
        <p className={styles.title}>{title}</p>
        
        {children && (
          <div className={styles.bioContainer}>
            {children}
          </div>
        )}
        
        <div className={styles.socialLinks}>
          {github && (
            <a 
              href={`https://github.com/${github}`} 
              target="_blank" 
              rel="noopener noreferrer"
              className={styles.socialLink}
              title={`GitHub: ${github}`}
            >
              <i className="fa fa-github"></i>
            </a>
          )}
          
          {linkedin && (
            <a 
              href={`https://linkedin.com/in/${linkedin}`} 
              target="_blank" 
              rel="noopener noreferrer"
              className={styles.socialLink}
              title={`LinkedIn: ${linkedin}`}
            >
              <i className="fa fa-linkedin" />
            </a>
          )}
          
          {twitter && (
            <a 
              href={`https://twitter.com/${twitter}`} 
              target="_blank" 
              rel="noopener noreferrer"
              className={styles.socialLink}
              title={`Twitter: ${twitter}`}
            >
              <i className="fa fa-twitter" />
            </a>
          )}
          
          {mastodon && (
            <a 
              href={mastodon} 
              target="_blank" 
              rel="noopener noreferrer"
              className={styles.socialLink}
              title="Mastodon"
            >
              <i className="fa fa-rss" />
            </a>
          )}
          
          {bluesky && (
            <a 
              href={bluesky} 
              target="_blank" 
              rel="noopener noreferrer"
              className={styles.socialLink}
              title="Bluesky"
            >
              <i className="fa fa-rss" />
            </a>
          )}
        </div>
      </div>
    </div>
  );
};

export default AmbassadorCard; 