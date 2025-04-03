import React, { useState, useEffect } from "react";
import clsx from "clsx";

import styles from "./styles.module.css";

type SlideProps = {
  className?: string;
  style?: React.CSSProperties;
  children: React.ReactNode;
};

const Slide: React.FC<SlideProps> = ({ className, style, children }) => {
  return (
    <div className={className} style={style}>
      {children}
    </div>
  );
};

type Props = {
  children: React.ReactNode;
  className?: string;
  width?: number;
  speed?: number;
};

const Marquee: React.FC<Props> = ({ children, className, width, ...props }) => {
  const [isMoving, setIsMoving] = useState(false);
  const [position, setPosition] = useState(0);

  const speed = props.speed || 30;

  useEffect(() => {
    setIsMoving(true);
    setPosition((pos) => pos - 100);
    const interval = setInterval(() => {
      setIsMoving(false);
      setTimeout(() => {
        setIsMoving(false);
        setPosition(0);
        setTimeout(() => {
          setIsMoving(true);
          setPosition(-100);
        }, 20);
      }, 0);
    }, speed * 1000);
    return () => clearInterval(interval);
  }, [speed]);

  return (
    <div className={clsx(styles.container, className)}>
      <Slide
        className={clsx(styles.slide1, { [styles.moving]: isMoving })}
        style={{
          transform: `translateX(${position}%)`,
          width: width ? `${width}px` : undefined,
          transitionDuration: isMoving ? `${speed}s` : undefined,
        }}
      >
        {children}
      </Slide>
      <Slide
        className={clsx(styles.slide2, { [styles.moving]: isMoving })}
        style={{
          transform: `translateX(${position}%)`,
          width: width ? `${width}px` : undefined,
          transitionDuration: isMoving ? `${speed}s` : undefined,
        }}
      >
        {children}
      </Slide>
    </div>
  );
};

export default Marquee;
