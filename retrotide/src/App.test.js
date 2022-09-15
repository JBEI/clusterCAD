import React from 'react';
import { render, screen } from '@testing-library/react';
import App from './App';

it('renders navigation bar', () => {
  render(<App />);
  const navElement = screen.getByText("RetroTide");
  expect(navElement).toBeInTheDocument();
});

it('renders correct copyright date', () => {
  const now = new Date();
  const thisYear = now.getFullYear();
  render(<App />);
  const copyYear = screen.queryByText(thisYear + " The Regents of", {exact: false});
  expect(copyYear).toBeInTheDocument();
});

it('fails to render wrong date', () => {
  const thisYear = 2929;
  render(<App />);
  const copyYear = screen.queryByText(thisYear + " The Regents of", {exact: false});
  expect(copyYear).not.toBeInTheDocument();
});