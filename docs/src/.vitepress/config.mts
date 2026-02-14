import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import { mathjaxPlugin } from './mathjax-plugin'
import footnote from "markdown-it-footnote";
import path from 'path'
import type { ShikiTransformer } from "shiki"

const mathjax = mathjaxPlugin()

type PromptKind = "julia" | "pkg" | "help" | "shell" | null

function juliaReplTransformer(): ShikiTransformer {
  let promptInfoByLine: Array<{ len: number; kind: PromptKind }> = []
  let isJuliaBlock = false

  function classify(line: string): { len: number; kind: PromptKind } {
    const rules: Array<{ kind: PromptKind; re: RegExp }> = [
      { kind: "julia", re: /^(julia>)(\s*)/ },
      { kind: "pkg", re: /^(\([^)]*\)\s*)?pkg>(\s*)/ },  // handles (@v1.9) pkg>
      { kind: "help", re: /^(help\?>)(\s*)/ },
      { kind: "shell", re: /^(shell>)(\s*)/ },
    ]

    for (const r of rules) {
      const m = line.match(r.re)
      if (m) return { len: m[0].length, kind: r.kind }
    }

    return { len: 0, kind: null }
  }

  return {
    name: "julia-repl-prompts",

    preprocess(code, options) {
      isJuliaBlock = options.lang === "julia"
    },

    tokens(tokens) {
      if (!isJuliaBlock) {
        promptInfoByLine = []
        return
      }
      
      promptInfoByLine = tokens.map((lineTokens) => {
        const line = lineTokens.map((t) => t.content).join("")
        return classify(line)
      })
    },

    span(node, line, col) {
      if (!isJuliaBlock) return
      
      const info = promptInfoByLine[line - 1]
      if (!info || !info.kind || info.len <= 0) return

      if (col < info.len) {
        this.addClassToHast(node, "repl-prompt")
        this.addClassToHast(node, `repl-prompt-${info.kind}`)
      }
    },
  }
}

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  lastUpdated: true,
  cleanUrls: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
  head: [
    ['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }],
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  
  markdown: {
    config(md) {
      md.use(tabsMarkdownPlugin);
      md.use(footnote);
      mathjax.markdownConfig(md);
    },
    codeTransformers: [juliaReplTransformer()],
    theme: {
      light: "github-light",
      dark: "github-dark"
    },
  },
  vite: {
    plugins: [
      mathjax.vitePlugin,
    ],
    define: {
      __DEPLOY_ABSPATH__: JSON.stringify('REPLACE_ME_DOCUMENTER_VITEPRESS_DEPLOY_ABSPATH'),
    },
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    optimizeDeps: {
      exclude: [ 
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ], 
    }, 
    ssr: { 
      noExternal: [ 
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ], 
    },
  },
  themeConfig: {
    outline: 'deep',
    logo: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    editLink: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
    socialLinks: [
      { icon: "slack", link: "https://julialang.org/slack/" },
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
