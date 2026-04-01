const { Document, Packer, Paragraph, TextRun, LevelFormat, AlignmentType } = require("docx");
const fs = require("fs");

const doc = new Document({
  styles: {
    default: { document: { run: { font: "Arial", size: 24 } } },
  },
  numbering: {
    config: [{
      reference: "bullets",
      levels: [{
        level: 0,
        format: LevelFormat.BULLET,
        text: "\u2022",
        alignment: AlignmentType.LEFT,
        style: { paragraph: { indent: { left: 720, hanging: 360 } } },
      }],
    }],
  },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },
        margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },
      },
    },
    children: [
      new Paragraph({
        spacing: { after: 300 },
        children: [new TextRun({ text: "Highlights", bold: true, size: 28 })],
      }),
      new Paragraph({
        numbering: { reference: "bullets", level: 0 },
        spacing: { after: 200 },
        children: [new TextRun({ text: "CFD-VisQA: first benchmark for evaluating VLM capability in CFD flow field plausibility assessment", size: 24 })],
      }),
      new Paragraph({
        numbering: { reference: "bullets", level: 0 },
        spacing: { after: 200 },
        children: [new TextRun({ text: "API-isolated protocol with base64-encoded images prevents filesystem contamination (99.6% vs 88.9%)", size: 24 })],
      }),
      new Paragraph({
        numbering: { reference: "bullets", level: 0 },
        spacing: { after: 200 },
        children: [new TextRun({ text: "Rank reversal: Claude leads with setup text (88.9%) but drops to worst without it (33.3%)", size: 24 })],
      }),
      new Paragraph({
        numbering: { reference: "bullets", level: 0 },
        spacing: { after: 200 },
        children: [new TextRun({ text: "Expert stability (\u22127.1 pp) reveals gestalt assessment grounded in physical intuition", size: 24 })],
      }),
      new Paragraph({
        numbering: { reference: "bullets", level: 0 },
        spacing: { after: 200 },
        children: [new TextRun({ text: "Setup-image cross-referencing, not visual physics understanding, drives VLM performance", size: 24 })],
      }),
    ],
  }],
});

const outPath = process.argv[2] || "highlights.docx";
Packer.toBuffer(doc).then(buffer => {
  fs.writeFileSync(outPath, buffer);
  console.log("Saved:", outPath);
});
