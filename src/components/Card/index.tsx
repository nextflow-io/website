// website/src/components/Card/index.tsx
import React from 'react';
import clsx from 'clsx';
import style from './styles.module.css';
import Button from '@components/Button';
interface CardProps {
    icon?: React.ReactNode; 
    title: string;          
    text: string;          
    link?: string;         
    isHoverable?: boolean; 
    button?: string;
}

const Card: React.FC<CardProps> = ({ icon, title, text, link, isHoverable, button }) => {
    const cardContent = (
        <div className={style.cardContent}>
            <div className={style.cardContentInner }>
            {icon && <div className={style.cardIcon}>{icon}</div>}
            <h3 className={style.cardTitle}>{title}</h3>
            <p className={style.cardText}>{text}</p>    
            </div>

            {button && <Button url={button} variant="link">Get Started</Button>}
        </div>
    );


    return link ? (
        <a href={link} className={clsx(style.card, { 'hoverable': isHoverable })}>
            {cardContent}
        </a>
    ) : (
        <div className={clsx(style.card, { 'hoverable': isHoverable })}>
            {cardContent}
        </div>
    );
};

export default Card;