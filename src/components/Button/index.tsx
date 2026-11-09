import React from 'react';
import clsx from 'clsx';
import styles from './styles.module.css';
import arrow from './src/arrow.svg';
interface ButtonProps {
    url: string; 
    variant?: 'primary' | 'secondary' | 'link'; 
    state?: 'default' | 'disabled'; 
    children: React.ReactNode;
    className?: string;
    target?: string;
}

function Button({ url, variant = 'primary', state = 'default', children, className, target }: ButtonProps) {
    return (
        <a
            href={url}
            target={target}
            className={clsx('button-link', styles.button, 'w-fit', styles[variant], styles[state], className)}
        >
            {children}
            {variant === 'link' && (
                <span className={styles.arrow}>
                    <img src={arrow.src} alt="arrow" />
                </span>
            )}
        </a>
    );
}

export default Button;