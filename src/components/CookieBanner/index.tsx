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
  const [isHidden, setIsHidden] = useState(true);

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
    const expireDate = new Date(Date.now() + 1000 * 60 * 60 * 24 * 365); // 1 año
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

  useEffect(() => {
    if (cookies.preferencesSet) {
      setIsHidden(true);
      if (cookies.preferencesChoice === "all") {
        giveConsent();
      }
    } else {
      setIsHidden(false);
    }
  }, [cookies]);

  if (isHidden) return null;

  return (
    <div className="fixed bottom-4 left-1/2 transform -translate-x-1/2 bg-gray-900 text-white p-4 rounded-lg flex items-center gap-4 shadow-lg z-50">
      <img src={CookieIcon.src} alt="Cookie" className="h-8 w-8 hidden sm:block" />
      <p className="text-sm">
        This website uses cookies to offer you a better browsing experience.
      </p>
      <div className="flex gap-2">
        <button onClick={denyAll} className="text-white px-4 py-2 rounded-md hover:bg-gray-700 transition">
          Essential only
        </button>
        <button onClick={acceptAll} className="bg-green-500 text-white px-4 py-2 rounded-md hover:bg-green-600 transition">
          Accept all cookies
        </button>
      </div>
    </div>
  );
};

export default CookieBanner;
