import React, { useState, useCallback } from "react";
import clsx from "clsx";
import "./ambassador.css";

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
  title,}) => {
  const [showBio, setShowBio] = useState(false);

  const toggleBio = useCallback(() => {
    setShowBio((prev) => !prev);
  }, []);

  const handleKeyDown = useCallback(
    (e: React.KeyboardEvent) => {
      if (e.key === "Enter" || e.key === " ") {
        e.preventDefault();
        toggleBio();
      }
    },
    [toggleBio],
  );

  return (
    <div
      className={clsx(
        "bg-white rounded-xs border border-brand-opacity transition-all relative overflow-hidden h-full",
        className,
      )}
    >
      <div className="relative aspect-square overflow-hidden bg-brand-1000">
        <img src={`img/${img}`} className="absolute inset-0 w-full h-full object-cover" alt={`${name}`} />
      </div>

      <div className="flex flex-col ">
        <div className="p-4 pb-3">
          <h2 className="text-lg font-semibold text-brand-1100 my-0">{name}</h2>
          {title && <p className="text-sm text-brand-900 mt-1 mb-2">{title}</p>}

          <div className="flex items-center mb-0 min-h-[24px]">
            {children && (
              <button
                onClick={toggleBio}
                onKeyDown={handleKeyDown}
                type="button"
                tabIndex={0}
                className={clsx(
                  "text-nextflow-600 text-xs cursor-pointer flex items-center gap-1 hover:text-nextflow-800 focus:outline-none focus:ring-2 focus:ring-nextflow-600 focus:ring-opacity-50 rounded-sm p-1 -ml-1 transition-all duration-200",
                  showBio && "text-nextflow-800",
                )}
                aria-label={showBio ? "Hide biography" : "Show biography"}
                aria-expanded={showBio}
              >
                <span>{showBio ? "Hide biography" : "Show biography"}</span>
                <i className={`fa fa-chevron-${showBio ? "down" : "up"} transition-transform duration-300`}></i>
              </button>
            )}
          </div>

          {children && (
            <div
              className={clsx(
                "bio-bottomsheet absolute inset-x-0 bottom-0 z-20 bg-white rounded-t-lg shadow-lg transform transition-transform duration-500 ease-in-out",
                showBio ? "translate-y-0" : "translate-y-full",
              )}
              aria-hidden={!showBio}
            >
              <div className="p-5 pt-8 max-h-[99vh] overflow-y-auto ambassador-bio-content pb-15">
                {children}

                <button
                  onClick={toggleBio}
                  onKeyDown={handleKeyDown}
                  type="button"
                  tabIndex={0}
                  className={clsx(
                    "text-nextflow-600 text-xs cursor-pointer flex items-center gap-1 hover:text-nextflow-800 focus:outline-none focus:ring-2 focus:ring-nextflow-600 focus:ring-opacity-50 rounded-sm p-1 -ml-1 transition-all duration-200",
                    showBio && "text-nextflow-800",
                  )}
                  aria-label={showBio ? "Hide biography" : "Show biography"}
                  aria-expanded={showBio}
                >
                  <span>{showBio ? "Hide biography" : "Show biography"}</span>
                  <i className={`fa fa-chevron-${showBio ? "down" : "up"} transition-transform duration-300`}></i>
                </button>
              </div>
            </div>
          )}
        </div>

        <div className="px-4 pb-4 flex gap-4 items-center social-container bg-white">
          {github && (
            <a
              href={`https://github.com/${github}`}
              target="_blank"
              rel="noopener noreferrer"
              className="text-brand-900 hover:text-nextflow-600 transition-colors relative group focus:outline-none focus:ring-2 focus:ring-nextflow-600 focus:ring-opacity-50 rounded-sm"
              title={`GitHub: ${github}`}
            >
              <i className="fa fa-github text-lg"></i>
              <span className="absolute bottom-0 left-0 w-0 h-0.5 bg-nextflow-600 group-hover:w-full transition-all duration-300"></span>
            </a>
          )}

          {linkedin && (
            <a
              href={`https://linkedin.com/in/${linkedin}`}
              target="_blank"
              rel="noopener noreferrer"
              className="text-brand-900 hover:text-nextflow-600 transition-colors relative group focus:outline-none focus:ring-2 focus:ring-nextflow-600 focus:ring-opacity-50 rounded-sm"
              title={`LinkedIn: ${linkedin}`}
            >
              <i className="fa fa-linkedin text-lg" />
              <span className="absolute bottom-0 left-0 w-0 h-0.5 bg-nextflow-600 group-hover:w-full transition-all duration-300"></span>
            </a>
          )}

          {twitter && (
            <a
              href={`https://twitter.com/${twitter}`}
              target="_blank"
              rel="noopener noreferrer"
              className="text-brand-900 hover:text-nextflow-600 transition-colors relative group focus:outline-none focus:ring-2 focus:ring-nextflow-600 focus:ring-opacity-50 rounded-sm"
              title={`Twitter: ${twitter}`}
            >
              <i className="fa fa-twitter text-lg" />
              <span className="absolute bottom-0 left-0 w-0 h-0.5 bg-nextflow-600 group-hover:w-full transition-all duration-300"></span>
            </a>
          )}

          {mastodon && (
            <a
              href={mastodon}
              target="_blank"
              rel="noopener noreferrer"
              className="text-brand-900 hover:text-nextflow-600 transition-colors relative group focus:outline-none focus:ring-2 focus:ring-nextflow-600 focus:ring-opacity-50 rounded-sm"
              title="Mastodon"
            >
              <i className="fa fa-rss text-lg" />
              <span className="absolute bottom-0 left-0 w-0 h-0.5 bg-nextflow-600 group-hover:w-full transition-all duration-300"></span>
            </a>
          )}

          {bluesky && (
            <a
              href={bluesky}
              target="_blank"
              rel="noopener noreferrer"
              className="text-brand-900 hover:text-nextflow-600 transition-colors relative group focus:outline-none focus:ring-2 focus:ring-nextflow-600 focus:ring-opacity-50 rounded-sm"
              title="Bluesky"
            >
              <i className="fa fa-rss text-lg" />
              <span className="absolute bottom-0 left-0 w-0 h-0.5 bg-nextflow-600 group-hover:w-full transition-all duration-300"></span>
            </a>
          )}

          {country && (
            <div className="absolute right-4 z-10 w-4 h-4 rounded-full overflow-hidden">
              <img
                className="w-full h-full object-cover"
                src={`https://flagicons.lipis.dev/flags/4x3/${country}.svg`}
                width="30"
                height="30"
                alt={`${country} flag`}
              />
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default AmbassadorCard;
