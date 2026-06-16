#!/usr/bin/env node
/**
 * Generates /.well-known/agent-skills/index.json with SHA-256 digests
 * computed from the bytes served at each skill URL.
 */
import { createHash } from "node:crypto";
import { readFileSync, writeFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = dirname(fileURLToPath(import.meta.url));
const publicDir = join(__dirname, "..", "public");
const siteOrigin = "https://nextflow.io";

const localSkills = [
  {
    name: "nextflow-docs",
    type: "skill-md",
    description: "Navigate Nextflow documentation, examples, and training resources.",
    path: ".well-known/agent-skills/nextflow-docs/SKILL.md",
  },
];

function sha256OfBytes(bytes) {
  const hex = createHash("sha256").update(bytes).digest("hex");
  return `sha256:${hex}`;
}

function sha256OfFile(filePath) {
  const bytes = readFileSync(filePath);
  return sha256OfBytes(bytes);
}

const skills = localSkills.map((skill) => {
  const filePath = join(publicDir, skill.path);
  return {
    name: skill.name,
    type: skill.type,
    description: skill.description,
    url: `${siteOrigin}/${skill.path}`,
    digest: sha256OfFile(filePath),
  };
});

const index = {
  $schema: "https://schemas.agentskills.io/discovery/0.2.0/schema.json",
  skills,
};

const outputPath = join(publicDir, ".well-known", "agent-skills", "index.json");
writeFileSync(outputPath, `${JSON.stringify(index, null, 2)}\n`);
console.log(`Wrote ${outputPath} with ${skills.length} skill(s).`);
