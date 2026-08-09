import assert from "node:assert/strict";
import { access, readFile } from "node:fs/promises";

const readJson = async (path) => JSON.parse(await readFile(path, "utf8"));

const corePlugins = await readJson(".obsidian/core-plugins.json");
const communityPlugins = await readJson(".obsidian/community-plugins.json");
const periodicNotes = await readJson(".obsidian/plugins/periodic-notes/data.json");
const templater = await readJson(".obsidian/plugins/templater-obsidian/data.json");
const dailyTemplate = await readFile("Templates/Daily-TEMPLATE.md", "utf8");
const dailyTheme = await readFile(".obsidian/snippets/Daily Note Themes.css", "utf8");
const generalTweaks = await readFile(".obsidian/snippets/General Tweaks.css", "utf8");
const readme = await readFile("README.md", "utf8");
const onThisDayBase = await readFile("Bases/On This Day.base", "utf8");
const appearance = await readJson(".obsidian/appearance.json");

assert.equal(corePlugins["daily-notes"], false, "Core Daily Notes must be disabled");
for (const plugin of ["calendar", "periodic-notes", "templater-obsidian"]) {
  assert.ok(communityPlugins.includes(plugin), `${plugin} must be enabled`);
}

assert.deepEqual(periodicNotes.daily, {
  format: "YYYY-MM-DD",
  template: "Templates/Daily-TEMPLATE.md",
  folder: "Journal/Daily",
  enabled: true,
});

assert.equal(templater.trigger_on_file_creation_mode, "folder", "Daily notes must trigger Templater by folder");
assert.deepEqual(templater.folder_templates, [
  {
    folder: "Journal/Daily",
    template: "Templates/Daily-TEMPLATE.md",
  },
]);

for (const token of [
  "<%*",
  "moment(tp.file.title, \"YYYY-MM-DD\", true)",
  "hide-properties",
  "weekday-<% date.day() %>",
  "##### ❮ [[<% date.clone().subtract(1, \"day\").format(\"YYYY-MM-DD\") %>]]",
  "![[Bases/On This Day.base]]",
]) {
  assert.ok(dailyTemplate.includes(token), `Daily template is missing ${token}`);
}

const inspirationIndex = dailyTemplate.indexOf("## 工作灵感");
const progressIndex = dailyTemplate.indexOf("## 工作进展");
assert.ok(inspirationIndex >= 0, "Daily template is missing the inspiration section");
assert.ok(inspirationIndex < progressIndex, "Inspiration must appear before work progress");
assert.match(onThisDayBase, /this\.file\.name\.slice\(4\)/, "On This Day must match the current month and day");
assert.doesNotMatch(onThisDayBase, /file\.inFolder\(/, "On This Day must query the full vault");

for (let day = 0; day <= 6; day += 1) {
  assert.ok(dailyTheme.includes(`.weekday-${day}`), `Missing weekday-${day} theme`);
}
assert.match(dailyTheme, /\.daily :is\(h5, \.HyperMD-header\.HyperMD-header-5\)/, "Daily navigation must be centered");
assert.match(dailyTheme, /\.hide-properties/, "Daily note properties must be hidden");
assert.equal(appearance.cssTheme, "AnuPpuccin", "AnuPpuccin must be the active theme");
assert.ok(appearance.enabledCssSnippets.includes("General Tweaks"), "General Tweaks must be enabled");
assert.match(generalTweaks, /\.hide-properties/, "General Tweaks must support hidden properties");

await access(".obsidian/themes/AnuPpuccin/theme.css");
await access(".obsidian/themes/AnuPpuccin/manifest.json");

const activeConfig = JSON.stringify({ corePlugins, periodicNotes, templater });
assert.doesNotMatch(activeConfig, /dialy|[A-Za-z]:[\\/]/i, "Active config has a legacy path");

for (const folder of ["Daily", "Weekly", "Monthly", "Quarterly", "Yearly"]) {
  await access(`Journal/${folder}/.workplaner-keep`);
}

for (const phrase of [
  "Templater 自动触发",
  "依赖 Templater",
  "AnuPpuccin",
  "## 工作灵感",
  "Journal/Daily/YYYY-MM-DD.md",
]) {
  assert.ok(readme.includes(phrase), `README is missing: ${phrase}`);
}

console.log("Journal configuration validation passed.");
