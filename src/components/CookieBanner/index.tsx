import { useEffect, useState } from "react";
import { useCookies } from "react-cookie";
import clsx from "clsx";
import CookieIcon from "./src/cookie.svg"; 

declare global {
  interface Window {
    gtag?: (...args: any[]) => void;
    posthog?: { opt_in_capturing: () => void };
  }
}

const CookieBanner = () => {
  const [cookies, setCookie] = useCookies(["preferencesSet", "preferencesChoice"]);
  const [isHidden, setIsHidden] = useState<boolean>(() => {
    return cookies.preferencesSet ? true : false;
  });

  useEffect(() => {
    if (cookies.preferencesSet) {
      setIsHidden(true);
    } else {
      setIsHidden(false);
    }
  }, [cookies.preferencesSet]);

  const giveConsent = () => {
    if (window.gtag) {
      window.gtag("consent", "update", {
        ad_user_data: "granted",
        ad_personalization: "granted",
        ad_storage: "granted",
        analytics_storage: "granted",
      });
    }
    if (window.posthog) {
      window.posthog.opt_in_capturing();
    }
  };

  const acceptAll = () => {
    const expireDate = new Date(Date.now() + 1000 * 60 * 60 * 24 * 365);
    setCookie("preferencesSet", "true", { expires: expireDate, path: "/" });
    setCookie("preferencesChoice", "all", { expires: expireDate, path: "/" });
    giveConsent();
    setIsHidden(true);
  };

  const denyAll = () => {
    setCookie("preferencesSet", "true", { path: "/" });
    setCookie("preferencesChoice", "essential", { path: "/" });
    setIsHidden(true);
  };

  if (isHidden) return null;

  return (
    <div className="relative">
      <div
        className={clsx(
          "w-full p-4 lg:rounded-md",
          "bg-brand text-white fixed bottom-0 lg:left-1/2 lg:mb-6 lg:transform lg:-translate-x-1/2 lg:w-[70%] z-[99]",
        )}
      >
        <div className="flex items-center justify-between flex-wrap -my-2">
          <div className="flex items-center py-2 md:pr-4">
            {/* <div className="mr-4 hidden sm:block">
              <img src={CookieIcon.src} alt="Cookie" className="h-8 w-8"/>
            </div>      */}
            <div className="flex-auto text-xs ml-5">
              This website uses cookies to offer you a better browsing
              experience. <br className="hidden md:block" />
              Find out more on{" "}
              <a
                href="https://seqera.io/privacy-policy/#cookies"
                className="text-white underline"
              >
                how we use cookies
              </a>
              .
            </div>
          </div>
          <div className="py-2">
            <div className="flex flex-wrap -mx-1 -my-1">
              <div className="px-1 py-1">
                <button
                  className=" text-white px-4 py-2 rounded-md hover:bg-brand-1000 outline outline-white outline-1 transition"
                  onClick={denyAll}
                >
                  Essential only
                </button>
              </div>
              <div className="px-1 py-1">
                <button
                  onClick={acceptAll}
                  className="bg-nextflow-600 text-white px-4 py-2 rounded-md hover:bg-nextflow-700 transition"
                >
                  Accept all cookies
                </button>
              </div>
            </div>
          </div>
        </div>
      </div>
    </div>
  );
};

export default CookieBanner;
