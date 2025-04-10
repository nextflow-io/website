import React, { useState } from "react";
import "./_styles.css";
import Square from "./svg/square.svg?raw";
import SvgClipboard from "./Clipboard";

interface Note {
  text?: string;
  link?: {
    text: string;
    url: string;
  };
  parts?: Array<{
    text: string;
    link?: {
      text: string;
      url: string;
    };
  }>;
}

interface Step {
  title: string;
  subtitle?: string;
  main: string;
  code?: string;
  note?: Note;
  pill?: string;
}

const steps: Step[] = [
  {
    title: "Check prerequisites",
    subtitle: "Java 11 or later is required",
    main: "Make sure Java 11 or later is installed on your computer by using the command:",
    code: "java -version",
    note: {
      text: "Note: If you are having trouble installing or upgrading Java check out our documentation ",
      link: { text: "here", url: "/docs/java-installation" },
    },
    pill: "Step 1"
  },
  {
    title: "Set up",
    subtitle: "Dead easy to install",
    main: "Enter this command in your terminal:",
    code: "curl -s https://get.nextflow.io | bash",
    note: {
      parts: [
        {
          text: "Nextflow can also be installed from ",
          link: { text: "Bioconda", url: "https://bioconda.github.io/recipes/nextflow/README.html" },
        },
        {
          text: " or downloaded from ",
          link: { text: "GitHub", url: "https://github.com/nextflow-io/nextflow/releases" },
        },
        {
          text: ". See the ",
          link: { text: "installation docs", url: "https://www.nextflow.io/docs/latest/install.html" },
        },
        {
          text: " for more information.",
        },
      ],
    },
    pill: "Step 2"
  },
  {
    title: "Launch",
    subtitle: "Try a simple demo",
    main: "Run the classic Hello world by entering the following command:",
    code: "./nextflow run hello",
    note: {
      parts: [
        {
          text: "You'll want to make the binary executable and move it to a directory in your PATH. ",
        },
        {
          text: "See the ",
          link: { text: "installation docs", url: "https://www.nextflow.io/docs/latest/install.html" },
        },
        {
          text: ".",
        },
      ],
    },
    pill: "Step 3"
  },
];

const TerminalComponent: React.FC = () => {
  const [currentStep, setCurrentStep] = useState<number>(0);
  const [copySuccess, setCopySuccess] = useState<boolean>(false);

  const copyToClipboard = (text: string) => {
    navigator.clipboard.writeText(text).then(
      () => {
        setCopySuccess(true);
        setTimeout(() => setCopySuccess(false), 2000);
      },
      (err) => {
        console.error("Error al copiar: ", err);
      },
    );
  };

  const formatCodeAsTag = (code: string) => {
    if (!code.includes("=") && !code.includes("<") && !code.includes(">")) {
      if (code.startsWith("curl") || code.startsWith("./nextflow") || code.startsWith("java")) {
        const parts = code.split(" ");
        let formattedCode = `<span style="color: #22863a">${parts[0]}</span> `;

        if (parts.length > 1) {
          formattedCode += parts
            .slice(1)
            .map((part) => {
              if (part.startsWith("-")) {
                return `<span style="color: #032f62">${part}</span>`;
              } else if (part.startsWith("http") || part.includes(".com") || part.includes(".io")) {
                return `<span style="color: #0366d6">${part}</span>`;
              }
              return part;
            })
            .join(" ");
        }

        return formattedCode;
      }
      return code;
    }

    let formatted = code
      .replace(/&lt;(\/?[a-zA-Z0-9_-]+)&gt;|<(\/?[a-zA-Z0-9_-]+)>|(\w+)=/g, (match, tag1, tag2, attr) => {
        if (tag1) return `&lt;<span style="color: #22863a">${tag1}</span>&gt;`;
        if (tag2) return `&lt;<span style="color: #22863a">${tag2}</span>&gt;`;
        if (attr) return `<span style="color: #005cc5">${attr}</span>=`;
        return match;
      })
      .replace(/"([^"]*)"/g, (match, value) => {
        return `"<span style="color: #032f62">${value}</span>"`;
      })
      .replace(/\b(\d+)\b/g, (match, num) => {
        return `<span style="color: #005cc5">${num}</span>`;
      });

    return formatted;
  };

  return (
    <div className="terminal-wrapper">
      <div>
        <div className="terminal-title">
          <div className="title-with-icon flex flex-col items-start gap-2">
            <div className="flex flex-row items-center gap-2">
              <span className="icon" dangerouslySetInnerHTML={{ __html: Square }} />

              <div className="flex flex-col ">
                <h2 className="my-0">{steps[currentStep].title}</h2>
                {steps[currentStep].subtitle && <p className="subtitle">{steps[currentStep].subtitle}</p>}
              </div>
            </div>
          </div>
        </div>
        <div className="terminal-content">
          <p className="main-text">{steps[currentStep].main}</p>

          {steps[currentStep].code && (
            <div className="code-block">
              <div className="code-content">
                <code dangerouslySetInnerHTML={{ __html: formatCodeAsTag(steps[currentStep].code || "") }}></code>
                <button
                  className="copy-button"
                  onClick={() => copyToClipboard(steps[currentStep].code || "")}
                  title="Copy to clipboard"
                >
                  <SvgClipboard />
                  <span className={`copy-tooltip ${copySuccess ? "visible" : ""}`}>Copied!</span>
                </button>
              </div>
            </div>
          )}
        </div>
        {steps[currentStep].note && (
          <p className="note-text">
            {steps[currentStep].note.parts ? (
              steps[currentStep].note.parts.map((part, i) => (
                <React.Fragment key={i}>
                  {part.text}
                  {part.link && (
                    <a href={part.link.url} className="link">
                      {part.link.text}
                    </a>
                  )}
                </React.Fragment>
              ))
            ) : (
              <>
                {steps[currentStep].note.text}{" "}
                {steps[currentStep].note.link && (
                  <a href={steps[currentStep].note.link.url} className="link">
                    {steps[currentStep].note.link.text}
                  </a>
                )}
              </>
            )}
          </p>
        )}
      </div>

      <div className="terminal-buttons" data-index={currentStep.toString()}>
        {steps.map((step, index) => (
          <div
            key={index}
            className="pill"
            data-active={(index === currentStep).toString()}
            onClick={() => setCurrentStep(index)}
          >
            {step.pill}
          </div>
        ))}
      </div>
    </div>
  );
};

export default TerminalComponent;
