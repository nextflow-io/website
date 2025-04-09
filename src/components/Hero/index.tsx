import clsx from "clsx";
import styles from "./styles.module.css";
import BackgroundGrid from "./src/BackgroundGrid";

const Hero = ({ children, large }: { children: React.ReactNode; large?: boolean }) => {
  return (
    <div className={clsx(styles.hero, "hero-section w-full px-0")}>
      <div className={styles.bgSvgContainer}>
        <BackgroundGrid />
      </div>

      <div
        className={styles.gridOverlay}
        style={{
          [styles.fade]: true,
        }}
      >
          <div
            className={clsx(styles.cell)}
            aria-hidden="true"
          />
      </div>

      {children}

      <slot />
    </div>
  );
};

export default Hero;
