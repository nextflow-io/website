import React, { useEffect, useRef } from 'react';

const BackgroundGrid: React.FC = () => {
  const sparklesRef = useRef<Array<SVGPathElement | null>>([]);
  
  useEffect(() => {
    const sparks = sparklesRef.current.filter((el): el is SVGPathElement => el !== null);
    
    const getRandomOrder = () => {
      const indices = Array.from({ length: sparks.length }, (_, i) => i);
      return indices.sort(() => Math.random() - 0.5);
    };
    
    const animationDuration = 5000;
    const delayBetweenSparks = 2000;
    const totalCycleDuration = sparks.length * delayBetweenSparks + animationDuration;
    
    const animateSpark = (sparkEl: SVGPathElement, delay: number) => {
      const classes = sparkEl.className.baseVal.split(' ');
      const sparkNumber = classes.find((cls: string) => cls.startsWith('sparkle-'))?.split('-')[1];
      
      let keyframeName = '';
      const sparkId = sparkEl.id;
      
      if (sparkId.includes('spark-0') || sparkId.includes('spark-7')) {
        keyframeName = 'moveVerticalDownUp';
      } else if (sparkId.includes('spark-1') || sparkId.includes('spark-2')) {
        keyframeName = 'moveVerticalUpDown';
      } else if (sparkId.includes('spark-3') || sparkId.includes('spark-4') || 
                 sparkId.includes('spark-5') || sparkId.includes('spark-6')) {
        keyframeName = 'moveHorizontalRightLeft';
      }
      
      sparkEl.style.opacity = '0';
      
      setTimeout(() => {
        const animation = sparkEl.animate(
          getKeyframes(keyframeName),
          {
            duration: animationDuration,
            easing: 'ease-in-out',
            fill: 'forwards'
          }
        );
        
        sparkEl.style.opacity = '0.8';
        
        animation.onfinish = () => {
          sparkEl.style.opacity = '0';
        };
      }, delay);
    };
    
    const getKeyframes = (animationType: string) => {
      switch (animationType) {
        case 'moveHorizontalLeftRight':
          return [
            { transform: 'translateX(-100vw)' },
            { transform: 'translateX(100vw)' }
          ];
        case 'moveHorizontalRightLeft':
          return [
            { transform: 'translateX(100vw)' },
            { transform: 'translateX(-100vw)' }
          ];
        case 'moveVerticalUpDown':
          return [
            { transform: 'translateY(-100vh)' },
            { transform: 'translateY(100vh)' }
          ];
        case 'moveVerticalDownUp':
          return [
            { transform: 'translateY(100vh)' },
            { transform: 'translateY(-100vh)' }
          ];
        default:
          return [];
      }
    };
    
    const startAnimationSequence = () => {
      sparks.forEach(spark => {
        spark.style.opacity = '0';
      });
      
      const randomOrder = getRandomOrder();
      
      randomOrder.forEach((sparkIndex, i) => {
        const delay = i * delayBetweenSparks;
        animateSpark(sparks[sparkIndex], delay);
      });
      
      setTimeout(startAnimationSequence, totalCycleDuration);
    };
    
    startAnimationSequence();
    
    return () => {
      sparks.forEach(spark => {
        const animations = spark.getAnimations();
        animations.forEach(anim => anim.cancel());
      });
    };
  }, []);
  
  const setSparkleRef = (index: number) => (el: SVGPathElement | null) => {
    sparklesRef.current[index] = el;
  };
  
  return (
    <svg width="100%" height="100%" viewBox="0 0 1472 778" fill="none"
        xmlns="http://www.w3.org/2000/svg" preserveAspectRatio="xMidYMid slice">
        <defs>
            <linearGradient id="paint0_linear_907_740" x1="172.5" y1="750" x2="172.5" y2="824" gradientUnits="userSpaceOnUse">
                <stop stopColor="#0CAE8E"/>
                <stop offset="0.96" stopColor="white" stopOpacity="0"/>
            </linearGradient>
            <linearGradient id="paint1_linear_907_740" x1="557.5" y1="35" x2="557.5" y2="-14" gradientUnits="userSpaceOnUse">
                <stop stopColor="#0CAE8E"/>
                <stop offset="0.96" stopColor="white" stopOpacity="0"/>
            </linearGradient>
            <linearGradient id="paint2_linear_907_740" x1="941.5" y1="16" x2="941.5" y2="-38" gradientUnits="userSpaceOnUse">
                <stop stopColor="#0CAE8E"/>
                <stop offset="0.96" stopColor="white" stopOpacity="0"/>
            </linearGradient>
            <linearGradient id="paint3_linear_907_740" x1="12" y1="108.836" x2="-60" y2="108.836" gradientUnits="userSpaceOnUse">
                <stop stopColor="#0CAE8E"/>
                <stop offset="0.96" stopColor="white" stopOpacity="0"/>
            </linearGradient>
            <linearGradient id="paint4_linear_907_740" x1="12" y1="311.5" x2="-60" y2="311.5" gradientUnits="userSpaceOnUse">
                <stop stopColor="#0CAE8E"/>
                <stop offset="0.96" stopColor="white" stopOpacity="0"/>
            </linearGradient>
            <linearGradient id="paint5_linear_907_740" x1="12" y1="498.5" x2="-60" y2="498.5" gradientUnits="userSpaceOnUse">
                <stop stopColor="#0CAE8E"/>
                <stop offset="0.96" stopColor="white" stopOpacity="0"/>
            </linearGradient>
            <linearGradient id="paint6_linear_907_740" x1="12" y1="682.5" x2="-60" y2="682.5" gradientUnits="userSpaceOnUse">
                <stop stopColor="#0CAE8E"/>
                <stop offset="0.96" stopColor="white" stopOpacity="0"/>
            </linearGradient>
            <linearGradient id="paint7_linear_907_740" x1="941.5" y1="749" x2="941.5" y2="823" gradientUnits="userSpaceOnUse">
                <stop stopColor="#0CAE8E"/>
                <stop offset="0.96" stopColor="white" stopOpacity="0"/>
            </linearGradient>
            <clipPath id="clip0_907_740">
                <rect y="0.335938" width="1472" height="777" rx="20" fill="white"/>
            </clipPath>
        </defs>
        <g id="BG_">
            <g clipPath="url(#clip0_907_740)">
                <g id="BG__2">
                    <line id="Line 8" x1="173.383" y1="0.335937" x2="173.383" y2="777.336" stroke="#160F26" strokeOpacity="0.15"/>
                    <line id="Line 17" x1="941.5" y1="35.002" x2="941.5" y2="777.002" stroke="#160F26" strokeOpacity="0.15"/>
                    <line id="Line 18" x1="1324.94" y1="0.335937" x2="1324.94" y2="777.336" stroke="#160F26" strokeOpacity="0.15"/>
                    <line id="Line 9" x1="-12.6836" y1="108.836" x2="1472" y2="108.836" stroke="#160F26" strokeOpacity="0.15"/>
                    <path id="Line 10" d="M0 218.336H1484.68" stroke="#160F26" strokeOpacity="0.15"/>
                    <path id="Line 13" d="M-12.6836 311.836H1472" stroke="#160F26" strokeOpacity="0.15"/>
                    <line id="Line 14" x1="-12.6836" y1="404.836" x2="1472" y2="404.836" stroke="#160F26" strokeOpacity="0.15"/>
                    <line id="Line 16" x1="-12.6836" y1="590.836" x2="1472" y2="590.836" stroke="#160F26" strokeOpacity="0.15"/>
                    <line id="Line 19" x1="-12.6836" y1="683.836" x2="1472" y2="683.836" stroke="#160F26" strokeOpacity="0.15"/>
                    <line id="Line 20" x1="-12.6836" y1="776.836" x2="1472" y2="776.836" stroke="#160F26" strokeOpacity="0.15"/>
                    <line id="Line 15" x1="-12.6836" y1="497.836" x2="1472" y2="497.836" stroke="#160F26" strokeOpacity="0.15"/>
                    <line id="Line 11" x1="557.5" y1="39.002" x2="557.5" y2="777.002" stroke="#160F26" strokeOpacity="0.15"/>
                </g>
                <path id="spark-0" ref={setSparkleRef(0)} className="sparkle sparkle-0" d="M173 750V824" stroke="url(#paint0_linear_907_740)" strokeOpacity="0.6" strokeWidth="2" strokeLinecap="round"/>
                <path id="spark-1" ref={setSparkleRef(1)} className="sparkle sparkle-1" d="M557 35V-14" stroke="url(#paint1_linear_907_740)" strokeOpacity="0.6" strokeWidth="2" strokeLinecap="round"/>
                <path id="spark-2" ref={setSparkleRef(2)} className="sparkle sparkle-2" d="M941 16L941 -38" stroke="url(#paint2_linear_907_740)" strokeOpacity="0.6" strokeWidth="2" strokeLinecap="round"/>
                <path id="spark-3" ref={setSparkleRef(3)} className="sparkle sparkle-3" d="M12 109.336L-60 109.336" stroke="url(#paint3_linear_907_740)" strokeOpacity="0.6" strokeWidth="2" strokeLinecap="round"/>
                <path id="spark-4" ref={setSparkleRef(4)} className="sparkle sparkle-4" d="M12 312L-60 312" stroke="url(#paint4_linear_907_740)" strokeOpacity="0.6" strokeWidth="2" strokeLinecap="round"/>
                <path id="spark-5" ref={setSparkleRef(5)} className="sparkle sparkle-5" d="M12 499L-60 499" stroke="url(#paint5_linear_907_740)" strokeOpacity="0.6" strokeWidth="2" strokeLinecap="round"/>
                <path id="spark-6" ref={setSparkleRef(6)} className="sparkle sparkle-6" d="M12 683L-60 683" stroke="url(#paint6_linear_907_740)" strokeOpacity="0.6" strokeWidth="2" strokeLinecap="round"/>
                <path id="spark-7" ref={setSparkleRef(7)} className="sparkle sparkle-7" d="M942 749V823" stroke="url(#paint7_linear_907_740)" strokeOpacity="0.6" strokeWidth="2" strokeLinecap="round"/>
            </g>
        </g>
    </svg>
  );
};

export default BackgroundGrid; 