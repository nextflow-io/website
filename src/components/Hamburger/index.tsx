import React from "react";
import styles from "./styles.module.css";

type Props = {
  toggleMenu: () => void;
  isOpen: boolean;
};

const Hamburger: React.FC<Props> = ({ toggleMenu, isOpen }) => {
  return (
    <button
      className={`${styles.hamburger} ${isOpen ? styles.open : ''}`}
      onClick={toggleMenu}
    >
      <span />
      <span />
      <span />
    </button>
  );
};

export default Hamburger; 