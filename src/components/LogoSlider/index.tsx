import React from "react";
import styles from "./styles.module.css";
import Marquee from "./Marquee";

import achilles from "./src/achilles.svg";
import astrazeneca from "./src/astrazeneca.svg";
import flagship from "./src/flagship.svg";


interface Props {}

const logos = [achilles, astrazeneca, flagship]; 

const LogoSlider: React.FC<Props> = () => {
  
  const customers = logos; // Usa el array de logos si no hay clientes

  return (
    <Marquee width={customers.length * 100} speed={80}>
      {customers.map((customer, index) => {
        if (!customer) return null; // Asegúrate de que el cliente no sea nulo
        return (
          <div key={index} className={styles.item}>
            <img
              src={customer.src}
              alt={`Logo ${index + 1}`} // Proporciona un alt descriptivo
              className="w-[70%] max-h-[40px]"
            />
          </div>
        );
      })}
    </Marquee>
  );
};

export default LogoSlider;