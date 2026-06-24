#!/usr/bin/env node
/**
 * Syncs skills from nextflow-io/agent-skills and generates
 * /.well-known/agent-skills/index.json.
 */
import { execSync } from "node:child_process";
import { cpSync, existsSync, mkdirSync, readdirSync, readFileSync, rmSync, statSync, writeFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const __dirname = dirname(fileURLToPath(import.meta.url));
const repoRoot = join(__dirname, "..");
const publicDir = join(repoRoot, "public");
const siteOrigin = "https://nextflow.io";
const agentSkillsRepo = "https://github.com/nextflow-io/agent-skills.git";
const cacheDir = join(repoRoot, ".cache", "agent-skills");
const skillsSourceDir = join(cacheDir, "skills");
const skillsDestDir = join(publicDir, ".well-known", "agent-skills");

function syncAgentSkillsRepo() {
  if (!existsSync(cacheDir)) {
    execSync(`git clone --depth 1 ${agentSkillsRepo} ${cacheDir}`, { stdio: "inherit" });
    return;
  }

  execSync("git fetch origin master && git reset --hard origin/master", {
    cwd: cacheDir,
    stdio: "inherit",
  });
}

function parseFrontmatter(content) {
  const match = content.match(/^---\r?\n([\s\S]*?)\r?\n---/);
  if (!match) {
    return {};
  }

  const fields = {};
  const lines = match[1].split("\n");
  let i = 0;

  while (i < lines.length) {
    const line = lines[i];
    const fieldMatch = line.match(/^([A-Za-z0-9_-]+):\s*(.*)$/);
    if (!fieldMatch) {
      i++;
      continue;
    }

    const key = fieldMatch[1];
    const value = fieldMatch[2];

    if (value === "|" || value === ">") {
      i++;
      const block = [];
      while (i < lines.length) {
        const blockLine = lines[i];
        if (/^[A-Za-z0-9_-]+:\s*/.test(blockLine)) {
          break;
        }
        if (blockLine.startsWith("  ")) {
          block.push(blockLine.slice(2));
        } else if (blockLine === "") {
          block.push("");
        } else {
          break;
        }
        i++;
      }
      fields[key] = block.join("\n").trim();
      continue;
    }

    fields[key] = value.trim();
    i++;
  }

  return fields;
}

function copySkills() {
  rmSync(skillsDestDir, { recursive: true, force: true });
  mkdirSync(skillsDestDir, { recursive: true });

  const skillDirs = readdirSync(skillsSourceDir).filter((entry) =>
    statSync(join(skillsSourceDir, entry)).isDirectory(),
  );

  const skills = [];

  for (const dirName of skillDirs.sort()) {
    const sourceFile = join(skillsSourceDir, dirName, "SKILL.md");
    if (!existsSync(sourceFile)) {
      continue;
    }

    const destDir = join(skillsDestDir, dirName);
    mkdirSync(destDir, { recursive: true });
    const destFile = join(destDir, "SKILL.md");
    cpSync(sourceFile, destFile);

    const content = readFileSync(destFile);
    const frontmatter = parseFrontmatter(content.toString("utf8"));
    const relativePath = `.well-known/agent-skills/${dirName}/SKILL.md`;

    skills.push({
      name: frontmatter.name ?? dirName,
      type: "skill-md",
      description: frontmatter.description ?? `Nextflow agent skill from nextflow-io/agent-skills (${dirName}).`,
      url: `${siteOrigin}/${relativePath}`,
    });
  }

  return skills;
}

syncAgentSkillsRepo();
const skills = copySkills();

const index = {
  $schema: "https://schemas.agentskills.io/discovery/0.2.0/schema.json",
  source: "https://github.com/nextflow-io/agent-skills",
  skills,
};

const outputPath = join(skillsDestDir, "index.json");
writeFileSync(outputPath, `${JSON.stringify(index, null, 2)}\n`);
console.log(`Wrote ${outputPath} with ${skills.length} skill(s) from nextflow-io/agent-skills.`);
