import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import { mathjaxPlugin } from './mathjax-plugin'
import { juliaReplTransformer } from './julia-repl-transformer'
import footnote from "markdown-it-footnote";
import path from 'path'

const mathjax = mathjaxPlugin()

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

const sidebar = {
  '/introduction/': navTemp.nav[0] ? [navTemp.nav[0]] : [],
  '/tutorial/': navTemp.nav[1] ? [navTemp.nav[1]] : [],
  '/how-to/': navTemp.nav[2] ? [navTemp.nav[2]] : [],
  '/explanation/': navTemp.nav[3] ? [navTemp.nav[3]] : [],
  '/reference/': navTemp.nav[4] ? [navTemp.nav[4]] : [],
}

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: 'REPLACE_ME_DOCUMENTER_VITEPRESS',// TODO: replace this in makedocs!
  title: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  description: 'REPLACE_ME_DOCUMENTER_VITEPRESS',
  lastUpdated: true,
  cleanUrls: true,
  ignoreDeadLinks: true,
  outDir: 'REPLACE_ME_DOCUMENTER_VITEPRESS', // This is required for MarkdownVitepress to work correctly...
  head: [
    ['link', { rel: 'icon', href: 'REPLACE_ME_DOCUMENTER_VITEPRESS_FAVICON' }],
    ['link', { rel: 'stylesheet', href: `${baseTemp.base}assets/center-images.css` }],
    ['link', { rel: 'stylesheet', href: `${baseTemp.base}assets/hide-documenter-example-output.css` }],
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
    siteTitle: false,
    logo: {
      light: '/logo.svg',
      dark: '/logo-dark.svg'
    },
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar,
    editLink: {
      text: 'Edit this page',
      pattern: ({ filePath }: { filePath: string }): string => {
        const base = 'https://github.com/JuliaHealth/KomaMRI.jl/edit/master'
        if (filePath.startsWith('tutorial/')) {
          const name = filePath.replace('tutorial/', '').replace('.md', '')
          return `${base}/examples/3.tutorials/lit-${name}.jl`
        }
        if (filePath.startsWith('tutorial-pluto/')) {
          const name = filePath.replace('tutorial-pluto/', '').replace('.md', '')
          return `${base}/examples/4.reproducible_notebooks/pluto-${name}.jl`
        }
        if (filePath.startsWith('explanation/')) {
          const name = filePath.replace('explanation/', '').replace('.md', '')
          const litPath = `${base}/docs/src/explanation/lit-${name}.jl`
          // TODO: this is a bit hacky, change to more robust solution
          // Only lit-generated pages (1-phantom, 2-motion) have a lit- source; others fall back to .md
          const litPages = ['1-phantom', '2-motion']
          if (litPages.includes(name)) return litPath
        }
        return `${base}/docs/src/${filePath}`
      }
    },
    socialLinks: [
      { icon: "slack", link: "https://julialang.org/slack/" },
    ],
    footer: {
      message:
        'Made with <a href="https://documenter.juliadocs.org/stable/" target="_blank"><strong>Documenter.jl</strong></a>, <a href="https://vitepress.dev" target="_blank"><strong>VitePress</strong></a>, and <a href="https://luxdl.github.io/DocumenterVitepress.jl/stable" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>Released under the MIT License. Powered by the <a href="https://www.julialang.org">Julia Programming Language</a>.<br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()} Koma Development Team.`,
    }
  }
})