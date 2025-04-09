import React, { useEffect, useState, useRef } from "react";
import "./menu.css";
import Hamburger from "../Hamburger";

const Menu = ({}) => {
  const [isVisible, setIsVisible] = useState(true);
  const [hasScrolled, setHasScrolled] = useState(false);
  const [isMenuOpen, setIsMenuOpen] = useState(false);
  const [isMobile, setIsMobile] = useState(false);
  const [activeDropdown, setActiveDropdown] = useState<number | null>(null);
  const lastScrollTop = useRef(0);
  const headerRef = useRef<HTMLDivElement>(null);
  const dropdownRefs = useRef<Array<HTMLLIElement | null>>([null, null, null]);

  useEffect(() => {
    const checkMobile = () => {
      setIsMobile(window.innerWidth <= 767);
    };
    checkMobile();
    window.addEventListener("resize", checkMobile);
    return () => {
      window.removeEventListener("resize", checkMobile);
    };
  }, []);

  useEffect(() => {
    setActiveDropdown(null);
  }, [isMobile]);

  useEffect(() => {
    if (isMobile) {
      const handleGlobalClick = (e: MouseEvent) => {
        if (activeDropdown !== null) {
          const activeDropdownElement = dropdownRefs.current[activeDropdown];
          if (activeDropdownElement && !activeDropdownElement.contains(e.target as Node)) {
            setActiveDropdown(null);
          }
        }
      };

      document.addEventListener('click', handleGlobalClick);
      return () => {
        document.removeEventListener('click', handleGlobalClick);
      };
    }
  }, [isMobile, activeDropdown]);

  useEffect(() => {
    if (!isMobile) {
      const handleDesktopClickOutside = (e: MouseEvent) => {
        if (activeDropdown !== null) {
          const activeDropdownElement = dropdownRefs.current[activeDropdown];
          if (activeDropdownElement && !activeDropdownElement.contains(e.target as Node)) {
            setActiveDropdown(null);
          }
        }
      };

      document.addEventListener('click', handleDesktopClickOutside);
      return () => {
        document.removeEventListener('click', handleDesktopClickOutside);
      };
    }
  }, [isMobile, activeDropdown]);

  useEffect(() => {
    const handleScroll = () => {
      const scrollTop = window.scrollY;

      if (scrollTop > 10 && !hasScrolled) {
        setHasScrolled(true);
      } else if (scrollTop <= 10 && hasScrolled) {
        setHasScrolled(false);
      }

      if (scrollTop > 200) {
        if (scrollTop > lastScrollTop.current) {
          setIsVisible(false);
        } else {
          setIsVisible(true);
        }
      } else {
        setIsVisible(true);
      }

      lastScrollTop.current = scrollTop <= 0 ? 0 : scrollTop;
    };

    if (!isMenuOpen) {
      window.addEventListener("scroll", handleScroll, { passive: true, capture: true });
      document.body.style.overflowY = "auto";
      handleScroll();
    } else {
      window.removeEventListener("scroll", handleScroll, { capture: true } as EventListenerOptions);
    }

    return () => {
      window.removeEventListener("scroll", handleScroll, { capture: true } as EventListenerOptions);
    };
  }, [hasScrolled, isMenuOpen]);

  const toggleMobileMenu = () => {
    setIsMenuOpen(!isMenuOpen);

    if (!isMenuOpen) {
      setIsVisible(true);
    }
  };

  const handleDropdownClick = (e: React.MouseEvent, index: number) => {
    e.preventDefault();
    e.stopPropagation();
    setActiveDropdown(activeDropdown === index ? null : index);
  };

  const handleMenuItemClick = (e: React.MouseEvent) => {
    e.stopPropagation();
  };

  const headerClasses = `navbar navbar-inverse navbar-fixed-top header-container
    ${!isVisible ? "header-hidden" : ""}
    ${isMenuOpen ? "menu-open" : ""}`;

  const handleLiClick = (e: React.MouseEvent, index: number) => {
    if (isMobile) {
      if ((e.target as HTMLElement).tagName === 'LI') {
        e.preventDefault();
        e.stopPropagation();
        setActiveDropdown(activeDropdown === index ? null : index);
      }
    }
  };

  const setDropdownRef = (index: number) => (el: HTMLLIElement | null) => {
    dropdownRefs.current[index] = el;
  };

  return (
    <div id="header" ref={headerRef} className={headerClasses} role="navigation">
      <div className="container ">
        <div className="navbar-header">
          {isMobile && (
            <div className="mobile-menu-toggle bottom-20">
              <Hamburger isOpen={isMenuOpen} toggleMenu={toggleMobileMenu} />
            </div>
          )}
          <a className="navbar-brand max-w-[155px] flex items-center" href="/index.html">
            <img src="/img/nextflow.svg" title="Nextflow Logo" className="mt-1" />
          </a>
        </div>
        <div className={`navbar-collapse ${isMenuOpen && isMobile ? "in" : "collapse"}`}>
          <ul className="nav navbar-nav py-2">
            <li className="show animated ">
              <a href="/docs/latest/index.html" className="text-black">
                Documentation
              </a>
            </li>

            <li className="show animated ">
              <a href="http://training.nextflow.io">
                Training
              </a>
            </li>

            <li className="show animated ">
              <a href="https://community.seqera.io/tag/nextflow" target="_blank" tabIndex={0}>
                Forum
              </a>
            </li>

            <li
              ref={setDropdownRef(0)}
              className={`dropdown show ${activeDropdown === 0 ? "open" : ""}`}
              onClick={(e) => handleLiClick(e, 0)}
            >
              <a
                href="#"
                className="dropdown-toggle"
                data-toggle="dropdown"
                aria-haspopup="true"
                aria-expanded={activeDropdown === 0 ? "true" : "false"}
                tabIndex={0}
                onClick={(e) => handleDropdownClick(e, 0)}
              >
                <span className="menu-label">Examples</span>
                <img src="/img/assets/angle-down.svg" alt="Expand" className="dropdown-icon inline-block" />
              </a>
              <ul className="dropdown-menu" role="menu" onClick={handleMenuItemClick}>
                <li>
                  <a href="/basic-pipeline.html" tabIndex={0}>
                    Basic pipeline
                  </a>
                </li>
                <li>
                  <a href="/mixing-scripting-languages.html" tabIndex={0}>
                    Mixing scripting languages
                  </a>
                </li>
                <li>
                  <a href="/blast-pipeline.html" tabIndex={0}>
                    BLAST pipeline
                  </a>
                </li>
                <li>
                  <a href="/rna-seq-pipeline.html" tabIndex={0}>
                    RNA-Seq pipeline
                  </a>
                </li>
                <li>
                  <a href="/machine-learning-pipeline.html" tabIndex={0}>
                    Machine Learning pipeline
                  </a>
                </li>
                <li>
                  <a href="https://github.com/nextflow-io/rnaseq-nf" target="_blank" tabIndex={0}>
                    Simple RNAseq pipeline
                  </a>
                </li>
                <li>
                  <a href="http://nextflow-io.github.io/patterns/index.html" tabIndex={0}>
                    Implementation patterns
                  </a>
                </li>
              </ul>
            </li>

            <li
              ref={setDropdownRef(1)}
              className={`dropdown show ${activeDropdown === 1 ? "open" : ""}`}
              onClick={(e) => handleLiClick(e, 1)}
            >
              <a
                href="#"
                className="dropdown-toggle"
                data-toggle="dropdown"
                aria-haspopup="true"
                aria-expanded={activeDropdown === 1 ? "true" : "false"}
                tabIndex={0}
                onClick={(e) => handleDropdownClick(e, 1)}
              >
                <span className="menu-label">Scientists</span>
                <img src="/img/assets/angle-down.svg" alt="Expand" className="dropdown-icon inline-block" />
              </a>
              <ul className="dropdown-menu" role="menu" onClick={handleMenuItemClick}>
                <li>
                  <a href="https://seqera.io/pipelines/" target="_blank" tabIndex={0}>
                    Pipelines
                  </a>
                </li>
                <li>
                  <a href="https://seqera.io/containers/" target="_blank" tabIndex={0}>
                    Containers
                  </a>
                </li>
                <li>
                  <a href="https://seqera.io/ask-ai/" target="_blank" tabIndex={0}>
                    Seqera AI
                  </a>
                </li>
              </ul>
            </li>



            <li
              ref={setDropdownRef(2)}
              className={`dropdown show ${activeDropdown === 2 ? "open" : ""}`}
              onClick={(e) => handleLiClick(e, 2)}
            >
              <a
                href="#"
                className="dropdown-toggle"
                data-toggle="dropdown"
                aria-haspopup="true"
                aria-expanded={activeDropdown === 2 ? "true" : "false"}
                tabIndex={0}
                onClick={(e) => handleDropdownClick(e, 2)}
              >
                <span className="menu-label">Resources</span>
                <img src="/img/assets/angle-down.svg" alt="Expand" className="dropdown-icon inline-block" />
              </a>
              <ul className="dropdown-menu" role="menu" onClick={handleMenuItemClick}>
                <li>
                  <a href="https://seqera.io/blog/tag-nextflow/" tabIndex={0}>
                    Blog
                  </a>
                </li>
                <li>
                  <a href="https://seqera.io/podcasts/" tabIndex={0}>
                    Podcast
                  </a>
                </li>
                <li>
                  <a href="https://community.seqera.io/tag/nextflow" tabIndex={0}>
                    Community forum
                  </a>
                </li>
                <li>
                  <a href="https://www.nextflow.io/slack-invite.html" tabIndex={0}>
                    Slack community chat
                  </a>
                </li>
                <li>
                  <a href="https://nf-co.re" tabIndex={0}>
                    nf-core pipelines
                  </a>
                </li>
                <li>
                  <a href="/ambassadors.html" tabIndex={0}>
                    Nextflow Ambassadors
                  </a>
                </li>
                <li>
                  <a href="/about-us.html" tabIndex={0}>
                    About Nextflow
                  </a>
                </li>
              </ul>
            </li>

            {isMobile && (
              <li className="navbar-right">
                <a href="https://github.com/nextflow-io/nextflow" title="GitHub Repository" tabIndex={0}>
                  <i className="fa fa-github hidden-xs" aria-hidden="true"></i>
                  <span className="visible-xs">GitHub repository</span>
                </a>
              </li>
            )}
          </ul>
          {!isMobile && (
          <ul className="navbar-right">
            <li className="navbar-right">
              <a href="https://github.com/nextflow-io/nextflow" title="GitHub Repository" tabIndex={0}>
                <i className="fa fa-github hidden-xs" aria-hidden="true"></i>
                <span className="visible-xs">GitHub repository</span>
              </a>
              </li>
            </ul>
          )}
        </div>
      </div>
    </div>
  );
};

export default Menu;
