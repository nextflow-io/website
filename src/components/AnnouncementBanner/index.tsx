import React, { useEffect, useRef } from "react";
import Button from "../Button";
import "./banner.css";

const FRONTIER_URL = "https://seqera.io/frontier/";

const AnnouncementBanner = () => {
  const bannerRef = useRef<HTMLDivElement>(null);

  // Keep the body offset in sync with the banner's real height (it wraps on small screens)
  useEffect(() => {
    const banner = bannerRef.current;
    if (!banner) return;
    const root = document.documentElement;
    const observer = new ResizeObserver(([entry]) => {
      root.style.setProperty("--banner-height", `${entry.target.getBoundingClientRect().height}px`);
    });
    observer.observe(banner);
    return () => {
      observer.disconnect();
      root.style.removeProperty("--banner-height");
    };
  }, []);

  return (
    <div ref={bannerRef} className="announcement-banner bg-brand" role="region" aria-label="Announcement">
      <div className="container announcement-banner__inner">
        <p className="announcement-banner__text">
          Our new product, <strong>Nextflow Frontier</strong>, is live, built on Nextflow by the team who created it.
        </p>
        <Button url={FRONTIER_URL} variant="link" className="announcement-banner__cta">
          Learn more
        </Button>
      </div>
    </div>
  );
};

export default AnnouncementBanner;
