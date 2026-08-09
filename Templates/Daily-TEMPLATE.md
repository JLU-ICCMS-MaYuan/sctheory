<%*
const date = moment(tp.file.title, "YYYY-MM-DD", true);
if (!date.isValid()) {
  throw new Error("Daily note filenames must use YYYY-MM-DD.");
}
const week = date.isoWeek().toString().padStart(2, "0");
const quarter = `Q${Math.floor(date.month() / 3) + 1}`;
-%>
---
date: "<% tp.file.title %>"
week: "[[Journal/Weekly/<% date.format("YYYY") %>-W<% week %>]]"
cssclasses:
  - hide-properties
  - daily
  - weekday-<% date.day() %>
---

## [[Journal/Yearly/<% date.format("YYYY") %>|<% date.format("YYYY") %>]] / [[Journal/Quarterly/<% date.format("YYYY") %>-<% quarter %>|<% quarter %>]] / [[Journal/Monthly/<% date.format("YYYY-MM") %>|<% date.format("MMMM") %>]] / [[Journal/Weekly/<% date.format("YYYY") %>-W<% week %>|Week <% date.isoWeek() %>]]

# DAILY NOTE

##### ❮ [[<% date.clone().subtract(1, "day").format("YYYY-MM-DD") %>]] | <% tp.file.title %> | [[<% date.clone().add(1, "day").format("YYYY-MM-DD") %>]] ❯

---

## 工作灵感

### 碎片想法

---

## 工作进展

### 目标

### 进展与产出

### 问题与决策

### 下一步

---

## 日记

### 碎片记录

### 情绪与反思

---

![[Bases/On This Day.base]]
